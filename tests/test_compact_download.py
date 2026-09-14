"""Workflow-only tests: no ROOT, credentials, EOS, or job submission."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

REPO=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("workflow",REPO/"runCondorAnalysis.py")
workflow=importlib.util.module_from_spec(spec);spec.loader.exec_module(workflow)


class CompactDownload(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.repo=Path(self.tmp.name);self.dest=self.repo/"rootfiles";self.dest.mkdir()
        campaign=self.repo/"condor/jobs/example_full";campaign.mkdir(parents=True)
        (campaign/"campaign.json").write_text(json.dumps({"jobs":[
            {"sample":s,"input_files":1} for s in ("data","mc","tt")]}))
        self.state={"workflow":"example","preset":"run2024i_legacy_tt","stage":"merged",
                    "full_campaign":"example_full","merged_directory":"root://example.invalid//eos/merged",
                    "compact_download":True,"merged_local_directory":str(self.dest)}
        self.path=self.repo/"state.json";self.calls=[]
        for s in ("DATA","MC","TT"):(self.dest/f"zjet_{s}.root").write_bytes(b"old")
        for p in (patch.object(workflow,"REPOSITORY",self.repo),patch.object(workflow,"confirm",return_value=True)):
            p.start();self.addCleanup(p.stop)

    def run_mock(self,cmd,**kw):
        self.calls.append(cmd)
        self.assertNotEqual(cmd[0],"xrdcp")
        if cmd[-1].startswith("writeCompactOutput.C"):
            target=Path(cmd[-1].split('"')[3]);target.write_bytes(b"compact")

    def test_three_compact_samples_without_full_download(self):
        with patch.object(workflow,"run",side_effect=self.run_mock):workflow.download_merged(self.path,self.state)
        self.assertEqual(self.state["stage"],"merged_downloaded")
        self.assertEqual(self.state["merged_local_content"],"histograms-only")
        self.assertEqual(len(self.calls),12)
        self.assertEqual(len(self.state["compact_sources"]),3)
        for s in ("DATA","MC","TT"):self.assertEqual((self.dest/f"zjet_{s}.root").read_bytes(),b"compact")
        self.assertFalse(list(self.dest.glob("*.part")))

    def test_failed_validation_preserves_original(self):
        def fail(cmd,**kw):
            self.run_mock(cmd)
            if cmd[-1].startswith("validateFlavorMatrix"):raise RuntimeError("bad compact")
        with patch.object(workflow,"run",side_effect=fail):
            with self.assertRaises(RuntimeError):workflow.download_merged(self.path,self.state)
        self.assertEqual(self.state["stage"],"merged")
        self.assertEqual((self.dest/"zjet_DATA.root").read_bytes(),b"old")
        self.assertFalse(list(self.dest.glob("*.part")))

    def test_reexport_completed_does_not_reset_checkpoint(self):
        self.state["stage"]="compatibility_written"
        with patch.object(workflow,"run",side_effect=self.run_mock):
            workflow.download_merged(self.path,self.state,advance_checkpoint=False)
        self.assertEqual(self.state["stage"],"compatibility_written")


if __name__=="__main__":unittest.main()
