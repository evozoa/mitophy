PY ?= python
PROFILE ?= ci

.PHONY: env pipeline pipeline-quick pipeline-deep site serve test clean

env:
	micromamba create -y -f environment.yml -n mitophy || conda env create -f environment.yml

pipeline:
	$(PY) -m mitophy all --profile $(PROFILE)

pipeline-quick:
	$(PY) -m mitophy all --profile quick --workdir work-quick

pipeline-deep:
	$(PY) -m mitophy all --profile deep

site:
	$(PY) -m mitophy site

serve: site
	$(PY) -m http.server -d _site 8000

test:
	$(PY) -m pytest -q

clean:
	rm -rf _site work-quick tests/work
