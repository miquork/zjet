"""Offline preparation only: never calls condor_submit or uses real credentials."""
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO/'scripts'))
import prepare_condor as prep

class TTPreparation(unittest.TestCase):
    @unittest.skipUnless((REPO/'data/Tagging/official_2024.txt').is_file() and
                        (REPO/'textfiles/generated/summer24_tt.txt').is_file(),
                        'Local private-derived WP configuration and generated lists required')
    def test_three_samples_and_private_config_transfer(self):
        with tempfile.TemporaryDirectory(prefix='zjet-tt-test-') as tmp:
            repo=Path(tmp)
            # Read-only source links; generated campaign lives only in the temp directory.
            for p in REPO.iterdir():
                if p.name not in ('condor','.git'): (repo/p.name).symlink_to(p,target_is_directory=p.is_dir())
            (repo/'condor').mkdir()
            (repo/'condor/run_zjet_job.sh').symlink_to(REPO/'condor/run_zjet_job.sh')
            proxy=Path('/afs/cern.ch/work/test/mock_proxy')
            original=Path.is_file
            def is_file(p):return True if p==proxy else original(p)
            argv=['prepare_condor.py','--campaign','test_tt',
                '--mc-list',str(REPO/'textfiles/generated/local_mc_test.txt'),
                '--data-list',str(REPO/'textfiles/generated/local_data_test.txt'),
                '--tt-list',str(REPO/'textfiles/generated/summer24_tt.txt'),
                '--max-mc-files','1','--max-data-files','1','--max-tt-files','1','--analysis-mode','legacy']
            with patch.object(prep,'REPOSITORY',repo),patch.object(sys,'argv',argv),\
                 patch.dict(os.environ,{'X509_USER_PROXY':str(proxy)}),patch.object(Path,'is_file',is_file):
                prep.main()
            campaign=repo/'condor/jobs/test_tt'
            metadata=json.loads((campaign/'campaign.json').read_text())
            self.assertEqual([j['sample'] for j in metadata['jobs']],['mc','data','tt'])
            self.assertEqual(metadata['tt_files'],1)
            self.assertEqual(metadata['analysis']['analysis_mode'],'legacy')
            self.assertEqual(metadata['jobs'][2]['output_file'],'zjet_TT_0000.root')
            submit=(campaign/'zjet.sub').read_text()
            self.assertIn(' legacy\n',submit)
            for source in ['ZJetTaggingControls.h','ZJetInputCounters.h','ZJetLegacyReplay.h','validateTaggingControls.C',
                           'data/Tagging/official_2024.txt']:
                self.assertIn(source,submit)
            self.assertNotIn('private/btagging.json',submit)

if __name__=='__main__':unittest.main()
