# Corpus query (run 2026-11-10)

PLOS Search API (Solr), paged with rows=1000 / start=0,1000,...:

```
GET https://api.plos.org/search
  q=subject:"Political science"
  fq=doc_type:full
  fq=journal_key:PLoSONE
  fq=publication_date:[2015-01-01T00:00:00Z TO *]
  fl=id,title_display,publication_date
  sort=id asc
  wt=json
```

numFound = 3,126 (the unrestricted subject query returns ~4,213; the date and journal
filters reduce it). All 3,126 DOIs were sorted lexicographically and 250 were sampled with
`random.Random(20261110).sample(...)` (see `fetch.py --n 250 --seed 20261110`).
The sample is listed in `candidates.tsv`.

JATS XML: `https://journals.plos.org/plosone/article/file?id=<DOI>&type=manuscript`
-> `raw/<paper_id>.xml` (4 concurrent requests, skip existing, one retry). 250/250 downloaded.
Note: the Solr index changes as PLOS publishes, so re-running the query later may yield a
different sample; `candidates.tsv` + `raw/` are the frozen record.
