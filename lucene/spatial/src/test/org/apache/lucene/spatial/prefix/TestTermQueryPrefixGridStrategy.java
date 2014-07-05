package org.apache.lucene.spatial.prefix;

/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import com.spatial4j.core.context.SpatialContext;
import com.spatial4j.core.shape.Shape;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.Field;
import org.apache.lucene.document.StoredField;
import org.apache.lucene.document.StringField;
import org.apache.lucene.spatial.SpatialTestCase;
import org.apache.lucene.spatial.prefix.tree.FlexPrefixTree2D;
import org.apache.lucene.spatial.prefix.tree.QuadPrefixTree;
import org.apache.lucene.spatial.prefix.tree.SpatialPrefixTree;
import org.apache.lucene.spatial.query.SpatialArgsParser;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;
import java.util.Arrays;


public class TestTermQueryPrefixGridStrategy extends SpatialTestCase {

  @Override
  @Before
  public void setUp() throws Exception {
    super.setUp();
    this.ctx = SpatialContext.GEO;
  }

  @Test
  public void testNGramPrefixGridLosAngeles_flex() throws  IOException {
    testNGramPrefixGridLosAngelesWithTrie(new FlexPrefixTree2D(ctx));
  }

  @Test
  public void testNGramPrefixGridLosAngeles_quad() throws  IOException {
    testNGramPrefixGridLosAngelesWithTrie(new QuadPrefixTree(ctx));
  }

  private void testNGramPrefixGridLosAngelesWithTrie(SpatialPrefixTree trie) throws IOException {
    TermQueryPrefixTreeStrategy prefixGridStrategy = new TermQueryPrefixTreeStrategy(trie, "geo");

    Shape point = ctx.makePoint(-118.243680, 34.052230);

    Document losAngeles = new Document();
    losAngeles.add(new StringField("name", "Los Angeles", Field.Store.YES));
    for (Field field : prefixGridStrategy.createIndexableFields(point)) {
      losAngeles.add(field);
    }
    losAngeles.add(new StoredField(prefixGridStrategy.getFieldName(), point.toString()));//just for diagnostics

    addDocumentsAndCommit(Arrays.asList(losAngeles));

    // This won't work with simple spatial context...
    SpatialArgsParser spatialArgsParser = new SpatialArgsParser();
    // TODO... use a non polygon query
//    SpatialArgs spatialArgs = spatialArgsParser.parse(
//        "Intersects(POLYGON((-127.00390625 39.8125,-112.765625 39.98828125,-111.53515625 31.375,-125.94921875 30.14453125,-127.00390625 39.8125)))",
//        new SimpleSpatialContext());

//    Query query = prefixGridStrategy.makeQuery(spatialArgs, fieldInfo);
//    SearchResults searchResults = executeQuery(query, 1);
//    assertEquals(1, searchResults.numFound);
  }
}
