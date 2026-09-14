"""Read-only reuse and provenance checks, without EOS or ROOT dependencies."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import merge_condor as merge


class MergeResume(unittest.TestCase):
    def test_provenance(self):
        fields=("campaign", "created_utc", "files_per_job", "mc_files", "data_files",
                "tt_files", "mc_jobs", "data_jobs", "tt_jobs", "source", "inputs", "analysis")
        expected={key:key for key in fields}
        actual=copy.deepcopy(expected)
        actual["merged_utc"]="old";actual["prepared_job_repairs"]=[]
        merge.compatible_provenance(actual,expected)
        for key in fields:
            wrong=copy.deepcopy(actual);wrong[key]="different"
            with self.assertRaisesRegex(RuntimeError,key): merge.compatible_provenance(wrong,expected)

    def test_local_reuse_never_changes_file(self):
        with tempfile.TemporaryDirectory() as d:
            scratch=Path(d);p=scratch/"zjet_MC.root";p.write_bytes(b"existing")
            with patch.object(merge,"validate_sample") as validate,patch.object(merge.subprocess,"run") as run:
                merge.reuse_existing(p,"mc",{}, {},scratch,False)
                validate.assert_called_once();run.assert_not_called()
            self.assertEqual(p.read_bytes(),b"existing")

    def test_failed_local_validation_leaves_file(self):
        with tempfile.TemporaryDirectory() as d:
            scratch=Path(d);p=scratch/"zjet_MC.root";p.write_bytes(b"existing")
            with patch.object(merge,"validate_sample",side_effect=RuntimeError("bad")):
                with self.assertRaises(RuntimeError): merge.reuse_existing(p,"mc",{}, {},scratch,False)
            self.assertEqual(p.read_bytes(),b"existing")

    def test_empty_file_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            scratch=Path(d);p=scratch/"zjet_MC.root";p.touch()
            with patch.object(merge,"validate_sample") as validate:
                with self.assertRaisesRegex(RuntimeError,"empty"):
                    merge.reuse_existing(p,"mc",{}, {},scratch,False)
                validate.assert_not_called()
            self.assertTrue(p.exists())

    def test_remote_reuse_downloads_only(self):
        with tempfile.TemporaryDirectory() as d:
            scratch=Path(d);url="root://example.invalid//eos/test/zjet_MC.root"
            def download(command,**kwargs):
                self.assertEqual(command[0:2],["xrdcp",url]);Path(command[2]).write_bytes(b"downloaded")
            with patch.object(merge.subprocess,"run",side_effect=download) as run,patch.object(merge,"validate_sample"):
                merge.reuse_existing(url,"mc",{}, {},scratch,True)
                run.assert_called_once()
            self.assertFalse((scratch/"existing_mc.root").exists())


if __name__=="__main__":unittest.main()
