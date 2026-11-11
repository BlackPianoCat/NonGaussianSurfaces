"""Packaging metadata: what users see on PyPI and after `pip install`."""
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
tomllib = pytest.importorskip('tomllib', reason='needs Python 3.11+') if sys.version_info >= (3, 11) else None
pytestmark = pytest.mark.skipif(not (ROOT / 'pyproject.toml').exists(), reason='needs a source checkout')


@pytest.fixture(scope='module')
def meta():
    if tomllib is None:
        pytest.skip('needs Python 3.11+')
    return tomllib.loads((ROOT / 'pyproject.toml').read_text())['project']


def test_names_and_entry_point(meta):
    from nongauss.cli import main
    assert meta['name'] == 'nongaussian-surfaces'
    assert meta['scripts']['nongauss'] == 'nongauss.cli:main' and callable(main)
    assert meta['readme'] == 'README.md' and (ROOT / 'README.md').exists()


def test_pypi_links_and_docs_url(meta):
    urls = meta['urls']
    assert {'Documentation', 'Source', 'Issues'} <= set(urls)
    cfg = (ROOT / 'mkdocs' / 'config.yml').read_text()
    assert f"site_url: {urls['Documentation']}" in cfg      # PyPI docs link = published docs site


def test_extras_and_install_names_in_docs(meta):
    assert {'interactive', 'docs', 'dev'} <= set(meta['optional-dependencies'])
    docs = (ROOT / 'mkdocs/docs/index.md').read_text() + (ROOT / 'README.md').read_text()
    assert f"pip install {meta['name']}" in docs
    assert f'pip install "{meta["name"]}[interactive]"' in docs


def test_publish_workflow_order():
    wf = ROOT / 'workflows' / 'publish.yml'
    wf = wf if wf.exists() else ROOT / '.github' / 'workflows' / 'publish.yml'
    text = wf.read_text()
    assert re.search(r'docs:\s*\n\s*needs: pypi', text)      # docs only after the PyPI upload
    assert re.search(r'pypi:\s*\n\s*needs: test', text)      # upload only after tests pass
    assert 'mkdocs gh-deploy' in text and 'mkdocs/config.yml' in text
