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

  private final int[] subCellsPerLevel;
  private final double[] levelW;
  private final double[] levelH;


  //The world bounds for the grid
  private final double xmin;
  private final double xmax;
  private final double ymin;
  private final double ymax;

  private final byte LEAF_BYTE = 0x01;
  private final byte START_CELL_NUMBER = 0x02;


  //The API will change- This is temporary for tests to pass
  public FlexPrefixTree2D(SpatialContext ctx, Rectangle bounds, int maxLevels) {
    this(ctx,bounds,maxLevels,null);
  }

  //Do we need maxlevels? SubcellsPerLevel should be enough
  //Can we provide a default cell division of 4
  private FlexPrefixTree2D(SpatialContext ctx, Rectangle bounds, int maxLevels,int[] subCellsPerLevel) {
    super(ctx,maxLevels);
    //Todo remove this
    if(subCellsPerLevel==null){
      subCellsPerLevel = new int[maxLevels+1];
      for(int i=0;i<=maxLevels;++i){
        subCellsPerLevel[i] = 2;
      }
    }

    this.subCellsPerLevel = subCellsPerLevel;

    this.xmin = bounds.getMinX();
    this.xmax = bounds.getMaxX();
    this.ymin = bounds.getMinY();
    this.ymax = bounds.getMaxY();

    levelW = new double[maxLevels+1];
    levelH = new double[maxLevels+1];

    double gridW = xmax - xmin;
    double gridH = ymax - ymin;

    //Now we will create a lookup for height and width for levels
    int division = this.subCellsPerLevel[0];
    levelW[0]= gridW;
    levelH[0]= gridH;

    //Compute the rest
    for (int i = 1; i < levelW.length; i++) {
      division = this.subCellsPerLevel[i];
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
  public int getLevelForDistance(double dist) {
    if (dist == 0)//short circuit
      return maxLevels;
    for (int i = 0; i < maxLevels-1; i++) {
      if(dist > levelW[i] && dist > levelH[i]) {
        return i;
      }
    }
    return maxLevels;
  }

  @Override
  public double getDistanceForLevel(int level) {
    if (level < 1 || level > getMaxLevels())
      throw new IllegalArgumentException("Level must be in 1 to maxLevels range");
    //get the grid width and height for that level
    double width = levelW[level];
    double height = levelH[level];
    //Use standard cartesian hypotenuse. For geospatial, this answer is larger
    // than the correct one but it's okay to over-estimate.
    return Math.sqrt(width * width + height * height);

  }

  private FlexCell[] newCellStack(int levels) {
    final FlexCell[] cellsByLevel = new FlexCell[levels + 1];
    final BytesRef term = new BytesRef(levels+1); //+1 For leaf and this byteRef will be shared within the stack
    for (int level = 0; level <= levels; level++) {
      cellsByLevel[level] = new FlexCell(term,cellsByLevel,level);
    }
    return cellsByLevel;
  }

  @Override
  public FlexCell getWorldCell() {
    return newCellStack(maxLevels)[0];
  }

  @Override
  public Cell readCell(BytesRef term, Cell scratch) {
    FlexCell cell = (FlexCell)scratch;
    if(scratch ==null){
      cell = getWorldCell();
    }
    //First get the length of the term
    int termLength = term.length;

    //We store at cell a cellstack len + leaf bytes
    //assert term.offset==0;
    termLength -= term.bytes[term.offset+term.length-1] == LEAF_BYTE?1:0;

    //Now from the cellstack obtain the correct numbered cell
    FlexCell cells[] = cell.getCellStack();
    cells[termLength].reuse(term);
    //assert cells[termLength].getLevel() == termLength;
    return cells[termLength];
  }

  /**
   * Cell here store the term without the leaf, so at any point we can remove
   * cellLevel and instead use term.length
   */
  private class FlexCell implements Cell {

    private final BytesRef term;

    protected SpatialRelation shapeRel;
    protected final FlexCell[] cellStack;
    private final int cellLevel;
    private final FlexPrefixTreeIterator cellIterator;

    private boolean isLeaf;
    private boolean isShapeSet= false;
    private final Rectangle gridRectangle;


    protected FlexCell(BytesRef term,FlexCell[] cells,int level){
      super();
      this.term = term;
      this.gridRectangle = ctx.makeRectangle(0,0,0,0);
      this.cellStack = cells;
      this.cellLevel = level;
      this.cellIterator = new FlexPrefixTreeIterator();
      readLeafAdjust();
    }

    private void readLeafAdjust() {
      isLeaf = (term.length > 0 && term.bytes[term.offset + term.length - 1] == LEAF_BYTE);
      if (isLeaf)
        term.length--;
    }

    protected FlexCell[] getCellStack(){
      return cellStack;
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
      this.isLeaf = true;
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
      term.length = cellLevel;
      if (result == null)
        return BytesRef.deepCopyOf(term);
      result.bytes = term.bytes;
      result.length = term.length;
      result.offset = term.offset;
      return result;
    }

    @Override
    public int getLevel() {
      return cellLevel;
    }

    @Override
    public CellIterator getNextLevelCells(Shape shapeFilter) {
      assert getLevel() < FlexPrefixTree2D.this.getMaxLevels();
      BytesRef source = getTokenBytesNoLeaf(null);
      int endCellNumber = (1<<subCellsPerLevel[getLevel()])+1;
      assert cellStack[getLevel()+1] != this: "Overwriting the parent";
      return cellIterator.initIter(source, cellStack[getLevel() + 1], shapeFilter, START_CELL_NUMBER, endCellNumber);
    }

    @Override
    public Shape getShape() {
      if(!isShapeSet) {
        isShapeSet =true;
        makeShape();
      }
      return gridRectangle;
    }

    private void makeShape() {
      BytesRef token = getTokenBytesNoLeaf(null);
      double xmin = FlexPrefixTree2D.this.xmin;
      double ymin = FlexPrefixTree2D.this.ymin;
      int col,row;
      boolean ancestorsSkippedRow=false;
      boolean ancestorsSkippedCol=false;
      //Todo avoid this loop by using the parent position and decoding only the bottom cell
      for (int i = 1; i <= token.length; i++) {
        int c = token.bytes[token.offset + i-1] -2;
        int division = subCellsPerLevel[i];
        col = (c / division);
        row = (c % division);
        xmin += levelW[i] * col;
        ymin += levelH[i] * row;
        if(row != division-1){
          ancestorsSkippedRow = true;
        }
        if(col != division-1){
          ancestorsSkippedCol = true;
        }
      }
      double xmax=xmin+levelW[token.length],ymax=ymin+levelH[token.length];
      if(token.length!=0) {
        //If parent has excluded the corner
        if (ancestorsSkippedRow) {
          ymax -= Math.ulp(ymax);
        }
        if (ancestorsSkippedCol) {
          xmax -= Math.ulp(xmax);
        }
      }
      gridRectangle.reset(xmin, xmax, ymin, ymax);
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

    private Cell reuse(BytesRef term) {
      this.isShapeSet = false;
      this.isLeaf = false;
      this.shapeRel = null;
      this.term.copyBytes(term);
      readLeafAdjust();
      return this;
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
    private int byte_pos;
    private int start;
    private int end;
    private FlexCell scratch;
    private Shape shapeFilter;

    //Inititalizes the Iterator, so that we can reuse the iterator
    protected CellIterator initIter(BytesRef source,FlexCell scratch,Shape shapeFilter,int start,int end){
      this.nextCell = null;
      this.thisCell = null;
      this.source = source;
      this.target = new BytesRef(source.length+1);
      this.target.copyBytes(source);
      this.byte_pos = this.target.length;
      this.target.length++;
      this.scratch = scratch;
      this.shapeFilter = shapeFilter;
      this.start = start;
      this.end = end;
      return this;
    }

    //Concatenates to the source BytesRef the given byte and places into te target
    private BytesRef concat(byte b) {
      assert target.offset == 0;
      target.bytes[byte_pos] = b;
      return target;
    }

    //Populates into scratch the next cell in z-order TODO Hilbert ordering
    private boolean hasNextCell(){
      if(start> end){
        nextCell=null;
        return false;
      }
      //We must call this as we want the cell to invalidate its ShapeCache
      scratch.reuse(concat((byte)start));
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
        assert nextCell.getLevel()==source.length+1: "The nextcell level is "+nextCell.getLevel()+"  source+1 "+(source.length+1);
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