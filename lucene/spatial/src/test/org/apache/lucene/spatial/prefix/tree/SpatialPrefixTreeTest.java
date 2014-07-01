package org.apache.lucene.spatial.prefix.tree;

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

import com.carrotsearch.randomizedtesting.annotations.Repeat;
import com.spatial4j.core.context.SpatialContext;
import com.spatial4j.core.context.SpatialContextFactory;
import com.spatial4j.core.shape.Point;
import com.spatial4j.core.shape.Rectangle;
import com.spatial4j.core.shape.Shape;
import com.spatial4j.core.shape.ShapeCollection;
import com.spatial4j.core.shape.SpatialRelation;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.Field;
import org.apache.lucene.document.Field.Store;
import org.apache.lucene.document.TextField;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.ScoreDoc;
import org.apache.lucene.search.TopDocs;
import org.apache.lucene.spatial.SpatialTestCase;
import org.apache.lucene.spatial.prefix.TermQueryPrefixTreeStrategy;
import org.apache.lucene.spatial.query.SpatialArgs;
import org.apache.lucene.spatial.query.SpatialOperation;
import org.junit.Before;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static com.carrotsearch.randomizedtesting.RandomizedTest.randomBoolean;
import static com.carrotsearch.randomizedtesting.RandomizedTest.randomIntBetween;

public class SpatialPrefixTreeTest extends SpatialTestCase {

  //TODO plug in others and test them
  private SpatialContext ctx;
  private SpatialPrefixTree trie;
  final int ITERATIONS = 10000;

  @Override
  @Before
  public void setUp() throws Exception {
    super.setUp();
    ctx = SpatialContext.GEO;
  }

  @Test
  public void testCellTraverse() {
    trie = new GeohashPrefixTree(ctx,4);

    Cell prevC = null;
    Cell c = trie.getWorldCell();
    assertEquals(0, c.getLevel());
    assertEquals(ctx.getWorldBounds(), c.getShape());
    while (c.getLevel() < trie.getMaxLevels()) {
      prevC = c;
      List<Cell> subCells = new ArrayList<>();
      CellIterator subCellsIter = c.getNextLevelCells(null);
      while (subCellsIter.hasNext()) {
        subCells.add(subCellsIter.next());
      }
      c = subCells.get(random().nextInt(subCells.size()-1));
      
      assertEquals(prevC.getLevel()+1,c.getLevel());
      Rectangle prevNShape = (Rectangle) prevC.getShape();
      Shape s = c.getShape();
      Rectangle sbox = s.getBoundingBox();
      assertTrue(prevNShape.getWidth() > sbox.getWidth());
      assertTrue(prevNShape.getHeight() > sbox.getHeight());
    }
  }
  /**
   * A PrefixTree pruning optimization gone bad, applicable when optimize=true.
   * See <a href="https://issues.apache.org/jira/browse/LUCENE-4770>LUCENE-4770</a>.
   */
  @Test
  public void testBadPrefixTreePrune() throws Exception {

    trie = new QuadPrefixTree(ctx, 12);
    TermQueryPrefixTreeStrategy strategy = new TermQueryPrefixTreeStrategy(trie, "geo");
    Document doc = new Document();
    doc.add(new TextField("id", "1", Store.YES));

    Shape area = ctx.makeRectangle(-122.82, -122.78, 48.54, 48.56);

    Field[] fields = strategy.createIndexableFields(area, 0.025);
    for (Field field : fields) {
      doc.add(field);
    }
    addDocument(doc);

    Point upperleft = ctx.makePoint(-122.88, 48.54);
    Point lowerright = ctx.makePoint(-122.82, 48.62);

    Query query = strategy.makeQuery(new SpatialArgs(SpatialOperation.Intersects, ctx.makeRectangle(upperleft, lowerright)));

    commit();

    TopDocs search = indexSearcher.search(query, 10);
    ScoreDoc[] scoreDocs = search.scoreDocs;
    for (ScoreDoc scoreDoc : scoreDocs) {
      System.out.println(indexSearcher.doc(scoreDoc.doc));
    }

    assertEquals(1, search.totalHits);
  }

  @Test
  @Repeat(iterations = ITERATIONS)
  public void testRandomCellRelationship(){

    int maxLevels = randomIntBetween(1,12);
    SpatialContextFactory ctxFactory = new SpatialContextFactory();
    ctxFactory.geo = false;
    ctxFactory.worldBounds = ctx.getWorldBounds();
    SpatialContext ctx = ctxFactory.newSpatialContext();
    assert ctx!= null;
    trie = new FlexPrefixTree2D(ctx,maxLevels);
    Rectangle WB = ctx.getWorldBounds();
    Point p = ctx.makePoint(randomIntBetween((int) WB.getMinX(), (int) WB.getMaxX()),randomIntBetween((int) WB.getMinY(), (int) WB.getMaxY()));
    //Get the world Cell
    Cell cell = trie.getWorldCell();
    CellIterator itr = cell.getNextLevelCells(null);

    if(maxLevels>1) {
      while (true) {
        //Get the cell that contains the point
        ArrayList<Shape> cells = new ArrayList<Shape>();
        Shape parent = cell.getShape();
        while (itr.hasNext()) {
          Cell c = itr.next();
          Rectangle s = (Rectangle) c.getShape();
          cells.add(ctx.makeRectangle(s.getMinX(), s.getMaxX(), s.getMinY(), s.getMaxY()));
        }
        Shape children = (new ShapeCollection<>(cells, ctx)).getBoundingBox();
        assertTrue(children.equals(parent));
        itr = cell.getNextLevelCells(p);
        cell = itr.next();
        if (cell.getLevel() >= (maxLevels - 2)) {
          break;
        }
        itr = cell.getNextLevelCells(null);
      }
    }
  }
/*
  @Test
  @Repeat(iterations = ITERATIONS)
  public void testRandomEdgeCellIntersectionsCount(){

    int maxLevels = randomIntBetween(1,12);
    SpatialContextFactory ctxFactory = new SpatialContextFactory();
    ctxFactory.geo = false;
    ctxFactory.worldBounds = ctx.getWorldBounds();
    SpatialContext ctx = ctxFactory.newSpatialContext();
    assert ctx!= null;
    trie = new FlexPrefixTree2D(ctx,maxLevels);
    Rectangle WB = ctx.getWorldBounds();
    double x = WB.getMaxX();
    double y = WB.getMinY();
    for(int i=1;i<maxLevels;++i){
      if(randomBoolean()){
        x /= 2; //Since FPT takes power of two
      }
      if(randomBoolean()){
        y /= 2;
      }
    }

    Point p= ctx.makePoint(x,y);
    Cell c = trie.getWorldCell();
    if(maxLevels>1) {
      checkNoOfMatches(p, c);
      //Check the centre of the grid TODO make it more points on the edge of the grid
      p = ctx.makePoint((WB.getMaxX() + WB.getMinX()) / 2, (WB.getMaxY() + WB.getMinY()) / 2);
      checkNoOfMatches(p, c);
    }
  }

  private void checkNoOfMatches(Point p,Cell c){
    CellIterator itr = c.getNextLevelCells(p);
    int count=0;
    while(itr.hasNext()){
      itr.next();
      ++count;
    }
    assertEquals(1,count);

  }*/
}