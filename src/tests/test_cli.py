import sys

import pytest

from nongauss.cli import main


@pytest.mark.parametrize('flag,ext', [([], 'png'), (['--interactive'], 'html')])
def test_cli_saves_files(tmp_path, monkeypatch, flag, ext):
    monkeypatch.setattr(sys, 'argv', ['nongauss', '--N', '48', '--out', str(tmp_path), *flag])
    main()
    assert len(list(tmp_path.glob(f'*.{ext}'))) >= 5
