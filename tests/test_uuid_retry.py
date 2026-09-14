"""Offline recovery-preparation tests: no Condor submission or credentials."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
import prepare_uuid_retry as retry


class RetryTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="zjet_uuid_retry_")
        self.addCleanup(self.tmp.cleanup)
        self.repo = Path(self.tmp.name)
        self.campaign = self.repo/"condor/jobs/test"
        self.campaign.mkdir(parents=True)
        self.logs = self.repo/"logs"; self.logs.mkdir()
        self.result = self.repo/"result.root"; self.result.touch()
        self.chunk = self.repo/"chunk.txt"; self.chunk.write_text("https://example.invalid/input.root\n")
        self.proxy = Path("/afs/test/protected_proxy")
        for name in ("zjet.C", "ZJetInputCounters.h", "worker.sh"):
            (self.repo/name).write_text("original\n")
        self.job = dict(sample="tt", chunk_id="0084", chunk_path="chunk.txt",
                        output_file="result.root", result_path=str(self.result), input_files=1)
        self.meta = dict(jobs=[self.job], storage=dict(log_directory=str(self.logs)),
                         source_files={name:dict(sha256=retry.digest(self.repo/name))
                                       for name in ("zjet.C", "ZJetInputCounters.h", "worker.sh")},
                         inputs=dict(tt_list=dict(selected_sha256=retry.digest(self.chunk))))
        (self.campaign/"campaign.json").write_text(json.dumps(self.meta))
        self.original = (f"initialdir = {self.repo}\nexecutable = worker.sh\n"
                         f"x509userproxy = {self.proxy}\ntransfer_input_files = zjet.C,ZJetInputCounters.h,$(chunk_path)\n"
                         f"output = {self.logs}/$(sample)_$(chunk_id).out\nlog = {self.logs}/condor.log\n"
                         f"queue sample,chunk_id,chunk_path,output_file,result_path from {self.campaign}/jobs.tsv\n")
        (self.campaign/"zjet.sub").write_text(self.original)
        (self.logs/"tt_0084.out").write_text("failed original log\n")
        (self.repo/"ZJetInputCounters.h").write_text("repaired\n")
        self.active = ""
        self.lifetime = "20000"
        def output(command, **kwargs):
            if command[0] == "condor_q": return self.active
            if command[-1] == "-timeleft": return self.lifetime
            if command[-1] == "-fqan": return "/cms/Role=NULL/Capability=NULL\n"
            raise AssertionError(command)
        is_file, stat = Path.is_file, Path.stat
        patches = [patch.object(retry, "__file__", str(self.repo/"scripts/prepare_uuid_retry.py")),
                   patch.object(retry.subprocess, "check_output", side_effect=output),
                   patch.object(Path, "is_file", lambda p: True if p == self.proxy else is_file(p)),
                   patch.object(Path, "stat", lambda p, **kw: SimpleNamespace(st_mode=0o100600)
                                if p == self.proxy else stat(p, **kw))]
        for p in patches: p.start(); self.addCleanup(p.stop)

    def prepare(self):
        return retry.prepare(str(self.campaign), "tt_0084", 1234, "uuidfix")

    def test_single_job_preserves_destinations_and_logs(self):
        directory = self.prepare()
        result = (directory/"zjet.sub").read_text()
        self.assertEqual(result.count("queue "), 1)
        self.assertIn("tt\t0084\tchunk.txt\tresult.root", result)
        self.assertIn(f"output = {self.logs}/$(sample)_$(chunk_id).out", result)
        self.assertIn(f"log = {directory}/condor.log", result)
        self.assertEqual((self.campaign/"zjet.sub").read_text(), self.original)
        self.assertEqual((directory/"original.out").read_text(), "failed original log\n")
        self.assertEqual((self.logs/"tt_0084.out").read_text(), "failed original log\n")
        record = json.loads((directory/"repair.json").read_text())
        self.assertEqual(set(record["changed_worker_sources"]), {"ZJetInputCounters.h"})
        with self.assertRaises(FileExistsError): self.prepare()

    def test_reject_nonempty(self):
        self.result.write_bytes(b"root")
        with self.assertRaisesRegex(ValueError, "nonempty"): self.prepare()

    def test_reject_active(self):
        self.active="1234 230\n"
        with self.assertRaisesRegex(ValueError, "queued"): self.prepare()

    def test_reject_expiring_proxy(self):
        self.lifetime="100"
        with self.assertRaisesRegex(ValueError, "Renew"): self.prepare()

    def test_reject_physics_change(self):
        (self.repo/"zjet.C").write_text("different physics\n")
        with self.assertRaisesRegex(ValueError, "physics"): self.prepare()

    def test_reject_changed_chunk(self):
        self.chunk.write_text("https://example.invalid/wrong.root\n")
        with self.assertRaisesRegex(ValueError, "selected-list hash"): self.prepare()

    def test_eos_schema_and_multiple_queues(self):
        original=self.original.replace(",result_path", "")
        result=retry.build_submit(original,self.job,self.campaign/"retry")
        self.assertIn("queue sample,chunk_id,chunk_path,output_file from (",result)
        with self.assertRaises(ValueError): retry.build_submit(original+"queue 1\n",self.job,self.campaign)


if __name__ == "__main__": unittest.main()
