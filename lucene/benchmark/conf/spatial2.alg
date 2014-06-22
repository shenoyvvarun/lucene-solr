#/**
# * Licensed to the Apache Software Foundation (ASF) under one or more
# * contributor license agreements.  See the NOTICE file distributed with
# * this work for additional information regarding copyright ownership.
# * The ASF licenses this file to You under the Apache License, Version 2.0
# * (the "License"); you may not use this file except in compliance with
# * the License.  You may obtain a copy of the License at
# *
# *     http://www.apache.org/licenses/LICENSE-2.0
# *
# * Unless required by applicable law or agreed to in writing, software
# * distributed under the License is distributed on an "AS IS" BASIS,
# * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# * See the License for the specific language governing permissions and
# * limitations under the License.
# */
# -------------------------------------------------------------------------------------
# Spatial search benchmark
#  In order to use this, you'll need to first run 'ant geonames-files'.
#  You may need more memory when running this: -Dtask.mem=1000M
#  For docs on what options are available, see the javadocs.

#  Benchmark existing geohash, and quad implementations
# -------------------------------------------------------------------------------------
### Spatial Context, Grid, Strategy config
#work.dir=work/geohash
#work.dir=work/quad
work.dir=work/flex
doc.maker=org.apache.lucene.benchmark.byTask.feeds.SpatialDocMaker
#spatial.prefixTree=geohash
#spatial.prefixTree=quad
spatial.prefixTree=flex
# 1m or better
spatial.maxDistErr = 0.000009
spatial.docPointsOnly=true
spatial.pruneLeafyBranches=false
#spatial.distErrPct=.25

### Source & Doc
content.source=org.apache.lucene.benchmark.byTask.feeds.LineDocSource
line.parser=org.apache.lucene.benchmark.byTask.feeds.GeonamesLineParser
docs.file=work/geonames/allCountries.txt
doc.tokenized=false

### Directory
directory=FSDirectory
compound=false
merge.factor=10
ram.flush.mb=64
concurrent.merge.scheduler.max.thread.count=2

### Query
query.maker=org.apache.lucene.benchmark.byTask.feeds.SpatialFileQueryMaker
query.file=work/geonames/allCountries.txt
query.file.line.parser=org.apache.lucene.benchmark.byTask.feeds.GeonamesLineParser
#query.file.maxQueries=1000

#  Next 3 props convert query points to circles with a random radius and then optionally bbox'es
query.spatial.radiusDegrees=0
query.spatial.radiusDegreesRandPlusMinus=3
query.spatial.bbox=false
query.spatial.prefixGridScanLevel=pgsl:-2:-3:-4:-5:-6:-7:-8:-9:-10:-11:-12
query.spatial.score=false

### Misc

log.step.AddDoc = 100000

{ "Populate"
  ResetSystemErase
  CreateIndex
  #1 million docs, to speed up add, do it parallely
  [{ "MAddDocs" AddDoc} : 250000] : 4
  ForceMerge(1)
  CommitIndex
  CloseIndex
} : 0

#set above round to 0 on subsequent runs if not changing indexing but experimenting with search

OpenReader
{"WarmJIT" Search > : 4000
CloseReader

{ "Rounds"
  ResetSystemSoft

  OpenReader
  Search
  {"Queries" Search > : 4000
  CloseReader

  NewRound
} : 10
RepSumByPrefRound Queries
