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

import com.spatial4j.core.context.SpatialContext;
import com.spatial4j.core.shape.Point;
import com.spatial4j.core.shape.Rectangle;
import com.spatial4j.core.shape.Shape;
import com.spatial4j.core.shape.SpatialRelation;
import org.apache.lucene.util.BytesRef;
import org.apache.lucene.util.StringHelper;

import java.util.ArrayList;
import java.util.List;

public class FlexPrefixTree2D extends SpatialPrefixTree{

  /**
   * Factory for creating {@link FlexPrefixTree2D} instances with useful defaults
   */
  public static class Factory extends SpatialPrefixTreeFactory {

    @Override
    protected int getLevelForDistance(double degrees) {
      FlexPrefixTree2D grid = new FlexPrefixTree2D(ctx, MAX_LEVELS_POSSIBLE);
      return grid.getLevelForDistance(degrees);
    }

    @Override
    protected SpatialPrefixTree newSPT() {
      return new FlexPrefixTree2D(ctx,
          maxLevels != null ? maxLevels : MAX_LEVELS_POSSIBLE);
    }
  }

  public static final int MAX_LEVELS_POSSIBLE = 50;//not really sure how big this should be

  public static final int DEFAULT_MAX_LEVELS = 12;

  public static final int START_CELL_NUMBER = 0x02;

  //The world bounds for the grid
  private final double xmin;
  private final double xmax;
  private final double ymin;
  private final double ymax;

  private final double gridW;
  public final double gridH;

  final double[] levelW; //width of each cell
  final double[] levelH;
  final int[] subcellsPerLevel;

  //The API will change- This is temporary for tests to pass
  public FlexPrefixTree2D(SpatialContext ctx, Rectangle bounds, int maxLevels) {
    this(ctx,bounds,maxLevels,null);
  }

  //Do we need maxlevels? SubcellsPerLevel should be enough
  //Can we provide a default cell division of 4
  private FlexPrefixTree2D(SpatialContext ctx, Rectangle bounds, int maxLevels,int []subcellsPerLevel) {
    super(ctx,maxLevels);
    //Todo remove this
    //If null is passed divide each level by 4
    if(subcellsPerLevel==null) {
      this.subcellsPerLevel = new int[maxLevels+1];
      for(int i=0;i<=maxLevels;++i){
        //TODO support odd power of 2
        // Everything works fine when its even power of 2
        //eg. 2^4= 16 cells can be configured as 4x4 cells
        //But, when its 2^3 = 8 cells, we can have either 2x4 or 4x2
        // Maybe alternate between 2x4 and 4x2 between levels
        this.subcellsPerLevel[i] = 2;
      }
      this.subcellsPerLevel[0] = 4;
    }else{
      this.subcellsPerLevel = subcellsPerLevel;
      assert this.subcellsPerLevel.length > 0 :"The subcells per Level must be greater than 0";
    }
    this.xmin = bounds.getMinX();
    this.xmax = bounds.getMaxX();
    this.ymin = bounds.getMinY();
    this.ymax = bounds.getMaxY();

    levelW = new double[maxLevels+1];
    levelH = new double[maxLevels+1];

    gridW = xmax - xmin;
    gridH = ymax - ymin;
    //Fill in the first element
    int division = this.subcellsPerLevel[0];// This is much much faster than division
    levelW[0]= gridW/division;
    levelH[0]=gridH/division;

    //Compute the rest
    for (int i = 1; i < levelW.length; i++) {
      division = this.subcellsPerLevel[i];
      levelW[i] = levelW[i - 1] / division;
      levelH[i] = levelH[i - 1] / division;
    }
  }

  //The API will change- This is temporary for tests to pass
  public FlexPrefixTree2D(SpatialContext ctx) {
    this(ctx, DEFAULT_MAX_LEVELS);
  }

  //The API will change- This is temporary for tests to pass
  public FlexPrefixTree2D(SpatialContext ctx, int maxLevels) {
    this(ctx, ctx.getWorldBounds(), maxLevels,null);
  }

  @Override
  public Cell getWorldCell() {
    return newCellStack(maxLevels)[0];
  }

  protected FlexCell[] newCellStack(int levels) {
    final FlexCell[] cellsByLevel = new FlexCell[levels + 1];
    final BytesRef term = new BytesRef(levels+1); //+1 For leaf and this byteRef will be shared within the stack
    for (int level = 0; level <= levels; level++) {
      cellsByLevel[level] = new FlexCell(cellsByLevel,term,level);
    }
    return cellsByLevel;
  }

  @Override
  public int getLevelForDistance(double dist) {
    if (dist == 0)//short circuit
      return maxLevels;
    for (int i = 0; i < maxLevels-1; i++) {
      //note: level[i] is actually a lookup for level i+1
      if(dist > levelW[i] && dist > levelH[i]) {
        return i+1;
      }
    }
    return maxLevels;
  }

  /**
   * Returns the cell containing point {@code p} at the specified {@code level}.
   */
  public Cell getCell(Point p, int level,FlexCell[] cellsByLevel) {
    //Since its a point there will be only a single cell match
    return build(xmin, ymin, 0, new BytesRef(maxLevels+1), ctx.makePoint(p.getX(),p.getY()), level,cellsByLevel);
  }

  //Given the bottom most corner it recursively checks over all cells at that particular level
  private Cell build(double x,double y,int level,BytesRef str,Shape shape,int maxLevel,FlexCell[] cellsByLevel) {
    assert str.length == level;
    double w = levelW[level];
    double h = levelH[level];
    int division = subcellsPerLevel[level];
    Cell match =null;
    // A z-ordering which is like an 'N' (a right rotated z)
    // Not a true Z order, but since we will be implementing Hilbert ordering later, now it will do
    // Writing this in a recursive manner will do the job
    cells: for(int i=0;i<division;++i){
      for(int j=0;j<division;++j){
        match = checkCellBehaviour((byte) ((i * division + j) + 2), x + (i * w), y + (j * h), level, str, shape, maxLevel, cellsByLevel);
        if(match != null){
          break cells;
        }
      }
    }
    return match;
  }

  private Cell checkCellBehaviour(byte c, double cx, double cy, int level, BytesRef str, Shape shape, int maxLevel, FlexCell[] cellsByLevel) {
    assert str.length == level;
    assert str.offset == 0;
    double w = levelW[level];
    double h = levelH[level];
    Cell match=null;
    int strlen = str.length;
    Rectangle rectangle = ctx.makeRectangle(cx,cx+w,cy,cy+h); // TODO reuse shape here
    SpatialRelation v = shape.relate(rectangle);
    if (SpatialRelation.CONTAINS == v) {
      str.bytes[str.length++] = c;//append
      match = cellsByLevel[level].reuse(str, v.transpose());
    } else if (SpatialRelation.DISJOINT == v) {
      // nothing
    } else { // SpatialRelation.WITHIN, SpatialRelation.INTERSECTS
      str.bytes[str.length++] = c;//append

      int nextLevel = level+1;
      if (nextLevel >= maxLevel) {
        match = cellsByLevel[level].reuse(str, v.transpose());
      } else {
        match = build(cx, cy, nextLevel,  str, shape, maxLevel,cellsByLevel);
      }
    }
    str.length = strlen;
    return match;
  }

  @Override
  public double getDistanceForLevel(int level) {
    if (level < 1 || level > getMaxLevels())
      throw new IllegalArgumentException("Level must be in 1 to maxLevels range");
    //TODO cache for each level
    FlexCell worldCell = (FlexCell)getWorldCell();
    Cell cell = getCell(ctx.getWorldBounds().getCenter(), level,worldCell.cellsByLevel);
    Rectangle bbox = cell.getShape().getBoundingBox();
    double width = bbox.getWidth();
    double height = bbox.getHeight();
    //Use standard cartesian hypotenuse. For geospatial, this answer is larger
    // than the correct one but it's okay to over-estimate.
    return Math.sqrt(width * width + height * height);
  }

  @Override
  public Cell readCell(BytesRef term, Cell scratch) {
    FlexCell cell = (FlexCell) scratch;
    if (cell == null)
      cell = (FlexCell) getWorldCell();
    cell.readCell(term);
    return cell;
  }

  public class FlexCell implements Cell {

    final private int cellLevel;
    final private BytesRef term;

    protected boolean isLeaf;

    /**
     * When set via getSubCells(filter), it is the relationship between this cell
     * and the given shape filter. Doesn't participate in shape equality.
     */
    protected SpatialRelation shapeRel;

    protected Shape shape;//cached
    private boolean isShapeCacheInvalidated=true;

    final FlexCell[] cellsByLevel;

    final FlexPrefixTreeIterator cellIterator;

    /** Warning: Refers to the same bytes (no copy). If {@link #setLeaf()} is subsequently called then it
     * may modify bytes. */
    FlexCell(FlexCell[] cellsByLevel,BytesRef term,int level) {
      super();
      this.cellIterator = new FlexPrefixTreeIterator();
      this.cellsByLevel = cellsByLevel;
      this.term = term;
      this.cellLevel = level;
      readLeafAdjust();
    }

    private Cell reuse(BytesRef term) {
      this.isShapeCacheInvalidated=true;
      this.isLeaf = false;
      this.shapeRel = null;
      this.term.copyBytes(term);
      readLeafAdjust();
      return this;
    }

    protected Cell reuse(BytesRef str, SpatialRelation shapeRel){
      this.reuse(str);
      this.shapeRel = shapeRel;
      return this;
    }

    protected FlexPrefixTree2D getGrid() { return FlexPrefixTree2D.this; }

    /**
     * Performant implementations are expected to implement this efficiently by
     * considering the current cell's boundary.
     * <p/>
     * Precondition: Never called when getLevel() == maxLevel.
     * Precondition: this.getShape().relate(p) != DISJOINT.
     */
    protected FlexCell getSubCell(Point p) {
      return (FlexCell) FlexPrefixTree2D.this.getCell(p, getLevel() + 1,cellsByLevel);
    }

    @Override
    public Shape getShape() {
      if(isShapeCacheInvalidated) {
        isShapeCacheInvalidated=false;
        makeShape();
      }
      return shape;
    }

    private void makeShape() {
      BytesRef token = getTokenBytesNoLeaf(null);
      double xmin = FlexPrefixTree2D.this.xmin;
      double ymin = FlexPrefixTree2D.this.ymin;
      int OverallColNo = 0;
      int OverallRowNo = 0;
      int BoundsRowNo = 0;

      for (int i = 0; i < token.length; i++) {
        int c = token.bytes[token.offset + i] -2;
        int division = subcellsPerLevel[i];
        int row = (c / division);
        int col = (c % division);
        xmin += levelW[i] * row;
        ymin += levelH[i] * col;
        OverallColNo += col;
        OverallRowNo += row;
        BoundsRowNo += (subcellsPerLevel[i]-1);
      }
      int len = token.length;
      double width, height;
      if (len > 0) {
        width = levelW[len-1];
        height = levelH[len-1];
      } else {
        width = gridW;
        height = gridH;
      }
      if(OverallColNo != BoundsRowNo){
        long temp = Double.doubleToLongBits(width);
        --temp;
        width = Double.longBitsToDouble(temp);
      }
      if(OverallRowNo != BoundsRowNo){
        long temp = Double.doubleToLongBits(height);
        --temp;
        height = Double.longBitsToDouble(temp);
      }

      if(shape == null) {
        shape = ctx.makeRectangle(xmin, xmin + width, ymin, ymin + height);
      }else{
        ((Rectangle)shape).reset(xmin, xmin + width, ymin, ymin + height);
      }

    }

    private static final byte LEAF_BYTE = 0x01;//NOTE: must sort before letters & numbers

    protected void readCell(BytesRef bytes) {
      shapeRel = null;
      this.isShapeCacheInvalidated = true;
      this.term.bytes = bytes.bytes;
      this.term.offset = bytes.offset;
      this.term.length = bytes.length;
      readLeafAdjust();
    }

    private void readLeafAdjust() {
      isLeaf = (term.length > 0 && term.bytes[term.offset + term.length - 1] == LEAF_BYTE);
      if (isLeaf)
        term.length--;
    }

    @Override
    public SpatialRelation getShapeRel() {
      return shapeRel;
    }

    @Override
    public void setShapeRel(SpatialRelation rel) {
      this.shapeRel = rel;
    }

    @Override
    public boolean isLeaf() {
      return isLeaf;
    }

    @Override
    public void setLeaf() {
      isLeaf = true;
    }

    @Override
    public BytesRef getTokenBytesWithLeaf(BytesRef result) {
      result = getTokenBytesNoLeaf(result);
      if (!isLeaf)
        return result;
      result.bytes[result.offset + result.length++] = LEAF_BYTE;
      return result;
    }

    @Override
    public BytesRef getTokenBytesNoLeaf(BytesRef result) {
      if (result == null)
        return BytesRef.deepCopyOf(term);
      result.bytes = term.bytes;
      result.length = term.length;
      result.offset = term.offset;
      return result;
    }

    @Override
    public int getLevel() {
      return term.length;
    }

    private int getCellLevel(){
      return cellLevel;
    }

    @Override
    public CellIterator getNextLevelCells(Shape shapeFilter) {
      assert getLevel() < getGrid().getMaxLevels();
      if (shapeFilter instanceof Point) {
        FlexCell cell = getSubCell((Point) shapeFilter);
        cell.shapeRel = SpatialRelation.CONTAINS;
        return new SingletonCellIterator(cell);
      } else {
        BytesRef source = getTokenBytesNoLeaf(null);
        int endCellNumber = (1<<subcellsPerLevel[getLevel()])+1;
        return cellIterator.initIter(source, cellsByLevel[getCellLevel() + 1], shapeFilter, START_CELL_NUMBER, endCellNumber);
      }
    }

    @Override
    public boolean isPrefixOf(Cell c) {
      //Note: this only works when each level uses a whole number of bytes.
      FlexCell cell = (FlexCell)c;
      boolean result = StringHelper.startsWith(cell.term, term);
      assert result == StringHelper.startsWith(c.getTokenBytesNoLeaf(null), getTokenBytesNoLeaf(null));
      return result;
    }

    @Override
    public int compareToNoLeaf(Cell fromCell) {
      FlexCell b = (FlexCell) fromCell;
      return b.term.compareTo(term);
    }

    @Override
    public boolean equals(Object obj) {
      //this method isn't "normally" called; just in asserts/tests
      if (obj instanceof Cell) {
        Cell cell = (Cell) obj;
        return getTokenBytesWithLeaf(null).equals(cell.getTokenBytesWithLeaf(null));
      } else {
        return false;
      }
    }

    @Override
    public String toString() {
      //this method isn't "normally" called; just in asserts/tests
      return getTokenBytesWithLeaf(null).utf8ToString();
    }

  }

  /**
   * An Iterator for FlexCells. This iterator reuses cells at a level and iterates over the siblings
   * initIter can be used to reuse the cell Iterator
   * TODO Use hilbert ordering
   */
  private class FlexPrefixTreeIterator extends CellIterator{

    private BytesRef source;
    private BytesRef target;
    private int start;
    private int end;
    private FlexCell scratch;
    private Shape shapeFilter;

    //Inititalizes the Iterator, so that we can reuse the iterator
    protected CellIterator initIter(BytesRef source,FlexCell scratch,Shape shapeFilter,int start,int end){
      this.nextCell = null;
      this.thisCell = null;
      this.source = source;
      this.target = scratch.term;
      this.scratch = scratch;
      this.shapeFilter = shapeFilter;
      this.start = start;
      this.end = end;
      return this;
    }

    //Concatenates to the source BytesRef the given byte and places into te target
    private BytesRef concat(BytesRef source, byte b, BytesRef target) {
      assert target.offset == 0;
      //target.bytes = new byte[source.length + 2];//+2 for new char + potential leaf
      target.length = 0;
      target.copyBytes(source);
      target.bytes[target.length++] = b;
      return target;
    }

    //Populates into scratch the next cell in z-order TODO Hilbert ordering
    private boolean hasNextCell(){
      if(start> end){
        return false;
      }
      scratch.reuse(concat(source, (byte)start, target), null);
      ++start;
      return true;
    }

    @Override
    public boolean hasNext() {
      thisCell = null;
      if (nextCell != null)//calling hasNext twice in a row
        return true;
      while (hasNextCell()) {
        nextCell = scratch;
        if (shapeFilter == null) {
          return true;
        } else {
          SpatialRelation rel = nextCell.getShape().relate(shapeFilter);
          if (rel.intersects()) {
            nextCell.setShapeRel(rel);
            if (rel == SpatialRelation.WITHIN)
              nextCell.setLeaf(); // Since the relation is a within no further decomposition will be required
            return true;
          }
        }
      }
      return false;
    }
  }
}
