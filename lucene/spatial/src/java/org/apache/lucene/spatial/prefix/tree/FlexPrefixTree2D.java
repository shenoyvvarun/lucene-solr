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

import java.util.Map;

import com.spatial4j.core.context.SpatialContext;
import com.spatial4j.core.shape.Point;
import com.spatial4j.core.shape.Rectangle;
import com.spatial4j.core.shape.Shape;
import com.spatial4j.core.shape.SpatialRelation;
import org.apache.lucene.util.BytesRef;
import org.apache.lucene.util.StringHelper;


public class FlexPrefixTree2D extends SpatialPrefixTree {

  /**
   * Factory for creating {@link FlexPrefixTree2D} instances with useful defaults
   */
  public static class Factory extends SpatialPrefixTreeFactory {

    public static final String LEVEL_PATTERN = "levelPattern";
    @Override
    protected void init(Map<String, String> args, SpatialContext ctx) {
      this.args = args;
      this.ctx = ctx;
      initNumberOfCellsPerLevel();
      initMaxLevels();
    }

    private void initNumberOfCellsPerLevel() {

    }

    @Override
    protected int getLevelForDistance(double degrees) {
      FlexPrefixTree2D grid = new FlexPrefixTree2D(ctx, MAX_LEVELS_POSSIBLE);
      return grid.getLevelForDistance(degrees);
    }

    @Override
    protected SpatialPrefixTree newSPT() {
      String levelPattern = args.get(LEVEL_PATTERN);
      if (levelPattern != null) {
        int[] numberOfSubCells;
        String[] levels = levelPattern.split(",");
        numberOfSubCells = new int[MAX_LEVELS_POSSIBLE];
        int i = 0;
        int repeatLevel=1;
        //Get the leading terms
        for (i = 0; i < levels.length; ++i) {
          if (levels[i].contains("*")) {
            levels[i] = levels[i].replace("*", "");
            repeatLevel = Integer.parseInt(levels[i]);
          }
          numberOfSubCells[i] = Integer.parseInt(levels[i]);
          maxLevels -= (numberOfSubCells[i]-1);
        }
        for (int j = i; j < MAX_LEVELS_POSSIBLE; ++j) {
          numberOfSubCells[j] = repeatLevel;
          maxLevels -= (numberOfSubCells[i]-1);
        }
        return new FlexPrefixTree2D(ctx, ctx.getWorldBounds(), maxLevels != null ? maxLevels : MAX_LEVELS_POSSIBLE, numberOfSubCells);
      } else {
        return new FlexPrefixTree2D(ctx, maxLevels != null ? maxLevels : MAX_LEVELS_POSSIBLE);
      }

    }
  }

  public static final int MAX_LEVELS_POSSIBLE = 30;//not really sure how big this should be

  public static final int DEFAULT_MAX_LEVELS = 12;

  private final int[] numberOfSubCellsAsExponentOfFour;
  private final int[] gridSizes;


  //The world bounds for the grid
  private final int xMin;
  private final int xMax;
  private final int yMin;
  private final int yMax;

  private final double intToDouble;
  private final double doubleToInt;
  private final double newOriginX;
  private final double newOriginY;
  private final Rectangle bounds;

  private final byte LEAF_BYTE = 0x01;
  private final byte START_CELL_NUMBER = 0x02;


  //The API will change- This is temporary for tests to pass
  public FlexPrefixTree2D(SpatialContext ctx, Rectangle bounds, int maxLevels) {
    this(ctx, bounds, maxLevels, null);
  }

  //Do we need maxlevels? SubcellsPerLevel should be enough
  //Can we provide a default cell division of 4
  public FlexPrefixTree2D(SpatialContext ctx, Rectangle bounds, int maxLevels, int[] numberOfSubCellsAsExponentOfTwo) {
    super(ctx, maxLevels);
    if (numberOfSubCellsAsExponentOfTwo == null) {
      numberOfSubCellsAsExponentOfTwo = new int[maxLevels];
      for (int i = 0; i < maxLevels; ++i) {
        numberOfSubCellsAsExponentOfTwo[i] = 1;
      }
    }
    this.bounds = bounds;
    this.numberOfSubCellsAsExponentOfFour = numberOfSubCellsAsExponentOfTwo;

    this.xMin = 0;
    this.yMin = 0;
    this.gridSizes = new int[maxLevels + 1];

    //Now we will create a lookup for height and width for levels
    int division;
    maxLevels = getMaxLevelsFromPowersOfFour(numberOfSubCellsAsExponentOfFour);
    assert maxLevels < 33;
    gridSizes[0] = (1 << maxLevels);

    //Compute the rest
    for (int i = 1; i < gridSizes.length; i++) {
      division = this.numberOfSubCellsAsExponentOfFour[i - 1];
      gridSizes[i] = (gridSizes[i - 1] >> division);
    }

    doubleToInt = Math.min((1 << maxLevels) / (bounds.getMaxY() - bounds.getMinY()), (1 << maxLevels) / (bounds.getMaxX() - bounds.getMinX()));

    intToDouble = Math.max((bounds.getMaxX() - bounds.getMinX()) / (1 << maxLevels), (bounds.getMaxY() - bounds.getMinY()) / (1 << maxLevels));

    newOriginX = bounds.getMinX();
    newOriginY = bounds.getMinY();
    xMax = (int) Math.ceil(bounds.getMaxX() * doubleToInt);
    yMax = (int) Math.ceil(bounds.getMaxX() * doubleToInt);
  }

  private int getMaxLevelsFromPowersOfFour(int[] numberOfCellsAsExponentOfFour) {
    int sum = 0;
    for (int i=0;i <maxLevels;++i) {
      sum += numberOfCellsAsExponentOfFour[i];
    }
    return sum;
  }

  //The API will change- This is temporary for tests to pass
  public FlexPrefixTree2D(SpatialContext ctx) {
    this(ctx, DEFAULT_MAX_LEVELS);
  }

  //The API will change- This is temporary for tests to pass
  public FlexPrefixTree2D(SpatialContext ctx, int maxLevels) {
    this(ctx, ctx.getWorldBounds(), maxLevels, null);
  }

  @Override
  public int getLevelForDistance(double dist) {
    if (dist == 0)//short circuit
      return maxLevels;
    for (int i = 1; i < maxLevels - 1; i++) {
      if (dist > (gridSizes[i]*intToDouble)) {
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
    double width = gridSizes[level] * intToDouble;
    double height = gridSizes[level] * intToDouble;
    //Use standard cartesian hypotenuse. For geospatial, this answer is larger
    // than the correct one but it's okay to over-estimate.
    return Math.sqrt(width * width + height * height);

  }

  @Override
  public FlexCell getWorldCell() {
    return new CellStack(maxLevels, xMin, yMin).cells[0];
  }

  @Override
  public Cell readCell(BytesRef term, Cell scratch) {
    FlexCell cell = (FlexCell) scratch;
    if (scratch == null) {
      cell = getWorldCell();
    }
    //First get the length of the term
    int termLength = term.length;

    //We store at cell a cellstack len + leaf bytes
    termLength -= term.bytes[term.offset + term.length - 1] == LEAF_BYTE ? 1 : 0;

    //Now from the cellstack obtain the correct numbered cell
    FlexCell cells[] = cell.getCellStack().cells;
    cells[termLength].reuseWithFreshTerm(term);
    cell.getCellStack().invalidate(0);
    return cells[termLength];
  }

  @Override
  public String toString() {
    StringBuilder subcells = new StringBuilder();
    for(int i=0;i<maxLevels;++i) {
      subcells.append(i);
      subcells.append("-");
    }
    return getClass().getSimpleName() + "(maxLevels:" + maxLevels + ",ctx:" + ctx +" ,NumberOfSubCells:"+subcells.toString()+")";
  }

  /**
   * Cell here store the term without the leaf, so at any point we can remove
   * cellLevel and instead use term.length
   */
  private class FlexCell implements Cell {

    protected final CellStack cellStack;
    private final int cellLevel;
    private final FlexPrefixTreeIterator cellIterator;
    private final Rectangle gridRectangle;
    private boolean isLeaf;
    private boolean isShapeSet = false;
    private SpatialRelation shapeRel;

    private int xMin;
    private int yMin;

    @Override
    public String toString() {
      return getShape() + "";
    }

    protected void setMinCornerCoordinates(int xmin, int ymin) {
      //Set the coordinates of bottom most corner
      this.xMin = xmin;
      this.yMin = ymin;
    }

    protected FlexCell(FlexCell cell, CellStack cells, int level) {
      super();
      this.gridRectangle = ctx.makeRectangle(0, 0, 0, 0);
      this.cellStack = cells;
      this.cellLevel = level;
      this.cellIterator = new FlexPrefixTreeIterator(cell, cellStack.term, cellLevel);
    }

    protected CellStack getCellStack() {
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
      cellStack.term.length = cellLevel;
      if (result == null)
        return BytesRef.deepCopyOf(cellStack.term);
      result.bytes = cellStack.term.bytes;
      result.length = cellStack.term.length;
      result.offset = cellStack.term.offset;
      return result;
    }

    @Override
    public int getLevel() {
      return cellLevel;
    }

    @Override
    public CellIterator getNextLevelCells(Shape shapeFilter) {
      assert getLevel() < FlexPrefixTree2D.this.getMaxLevels();
      return cellIterator.init(shapeFilter, START_CELL_NUMBER);
    }

    @Override
    public Shape getShape() {
      if (!isShapeSet) {
        isShapeSet = true;
        makeShape();
      }
      return gridRectangle;
    }

    private void makeShape() {
      BytesRef token = cellStack.term;
      cellStack.decode(cellLevel);//Decode the terms
      double xMax = ((xMin + gridSizes[token.length]) * intToDouble) + newOriginX;
      double yMax = ((yMin + gridSizes[token.length]) * intToDouble) + newOriginY;
      double xMin = (this.xMin * intToDouble) + newOriginX;
      double yMin = (this.yMin * intToDouble) + newOriginY;
      if (xMax < bounds.getMaxX()) {
        xMax = Math.nextAfter(xMax, Double.NEGATIVE_INFINITY);
      }

      if (yMax < bounds.getMaxY()) {
        yMax = Math.nextAfter(yMax, Double.NEGATIVE_INFINITY);
      }
      gridRectangle.reset(xMin, xMax, yMin, yMax);
    }

    @Override
    public boolean isPrefixOf(Cell c) {
      //Note: this only works when each level uses a whole number of bytes.
      FlexCell cell = (FlexCell) c;
      boolean result = StringHelper.startsWith(cell.cellStack.term, cellStack.term);
      assert result == StringHelper.startsWith(c.getTokenBytesNoLeaf(null), getTokenBytesNoLeaf(null));
      return result;
    }

    @Override
    public int compareToNoLeaf(Cell fromCell) {
      FlexCell b = (FlexCell) fromCell;
      return b.cellStack.term.compareTo(cellStack.term);
    }

    private Cell reuseWithFreshTerm(BytesRef term) {
      cellStack.invalidate(cellLevel); //Force re-decoding of the invalidated cell
      this.isShapeSet = false;
      this.isLeaf = false;
      this.shapeRel = null;
      this.cellStack.term.length = term.length;
      this.cellStack.term.bytes = term.bytes;
      this.cellStack.term.offset = term.offset;
      // Now the term placed here may have a leaf
      isLeaf = (cellStack.term.length > 0 && cellStack.term.bytes[cellStack.term.offset + cellStack.term.length - 1] == LEAF_BYTE);
      if (isLeaf)
        cellStack.term.length--;
      return this;
    }

    private Cell reuse() {
      cellStack.invalidate(cellLevel);
      this.isShapeSet = false;
      this.isLeaf = false;
      this.shapeRel = null;
      return this;
    }


    private SpatialRelation relateIntegerCoordinate(int int_min, int int_max, int ext_min, int ext_max) {
      if (ext_min > int_max || ext_max < int_min) {
        return SpatialRelation.DISJOINT;
      }

      if (ext_min >= int_min && ext_max <= int_max) {
        return SpatialRelation.CONTAINS;
      }

      if (ext_min <= int_min && ext_max >= int_max) {
        return SpatialRelation.WITHIN;
      }
      return SpatialRelation.INTERSECTS;
    }

    protected SpatialRelation relateIntegerRectangle() {
      int xMax = this.xMin + gridSizes[this.cellLevel];
      if (xMax < FlexPrefixTree2D.this.xMax) {
        xMax -= 1;
      }
      int ymax = this.yMin + gridSizes[this.cellLevel];
      if (ymax < FlexPrefixTree2D.this.yMax) {
        ymax -= 1;
      }
      SpatialRelation yIntersect = relateIntegerCoordinate(this.yMin, ymax, this.cellStack.shapeFilterYMin, this.cellStack.shapeFilterYMax);
      if (yIntersect == SpatialRelation.DISJOINT)
        return SpatialRelation.DISJOINT;
      SpatialRelation xIntersect = relateIntegerCoordinate(this.xMin, xMax, this.cellStack.shapeFilterXMin, this.cellStack.shapeFilterXMax);
      if (xIntersect == SpatialRelation.DISJOINT)
        return SpatialRelation.DISJOINT;
      if (xIntersect == yIntersect)//in agreement
        return xIntersect;
      if (this.cellStack.shapeFilterXMin == this.xMin && this.cellStack.shapeFilterXMax == xMax)
        return yIntersect;
      if (this.cellStack.shapeFilterYMin == this.yMin && this.cellStack.shapeFilterYMax == ymax)
        return xIntersect;
      return SpatialRelation.INTERSECTS;
    }

  }


  /**
   * An Iterator for FlexCells. This iterator reuses cells at a level and iterates over the siblings
   * initIter can be used to reuse the cell Iterator
   */
  private class FlexPrefixTreeIterator extends CellIterator {

    private final BytesRef term;
    private final int bytePos;
    private final int endCellNumber;
    private final FlexCell cell;
    private int nextCellNumber;
    private Shape shapeFilter;


    protected FlexPrefixTreeIterator(FlexCell cell, BytesRef sharedTerm, int level) {
      this.term = sharedTerm;
      this.cell = cell;
      if (level < maxLevels) {
        this.endCellNumber = (1 <<( numberOfSubCellsAsExponentOfFour[level] + numberOfSubCellsAsExponentOfFour[level])) + 1;
      } else {
        this.endCellNumber = 0;
      }
      this.bytePos = level;
    }

    //Inititalizes the Iterator, so that we can reuse the iterator
    protected CellIterator init(Shape shapeFilter, int start) {
      this.nextCell = null;
      this.thisCell = null;
      //Level 0 does not store a byte its byte pos is -1, but, in makeshape this is handled
      //Level 1 stores its byte at index 0
      //this.bytePos = this.scratch.cellLevel-1;
      this.shapeFilter = shapeFilter;
      this.cell.cellStack.findIntegerBoundingBox(shapeFilter);
      this.nextCellNumber = start;
      return this;
    }

    //Concatenates to the source BytesRef the given byte and places into te target
    private void changeTailByte(byte b) {
      term.bytes[term.offset + bytePos] = b;
    }

    @Override
    public boolean hasNext() {
      thisCell = null;
      if (nextCell != null)//calling hasNext twice in a row
        return true;
      while (levelHasUntraversedCell()) {
        SpatialRelation rel = null;
        nextCell = cell;
        if (shapeFilter == null) {
          return true;
        } else {
          FlexCell nextFlexCell = (FlexCell) nextCell;
          nextFlexCell.cellStack.decode(nextFlexCell.cellLevel);
          rel = getSpatialRelation(nextFlexCell);
          if (rel.intersects()) {
            nextCell.setShapeRel(rel);
            if (rel == SpatialRelation.WITHIN)
              nextCell.setLeaf(); // Since the relation is a within no further decomposition will be required
            if (rel == SpatialRelation.CONTAINS) {
              stopLevelIteration();
            }
            return true;
          }
        }
      }
      return false;
    }

    private SpatialRelation getSpatialRelation(FlexCell nextFlexCell) {
      SpatialRelation rel = null;
      if (!nextFlexCell.cellStack.shapeFilterBoundingBox.getCrossesDateLine()) {
        rel = nextFlexCell.relateIntegerRectangle();
        if (!(shapeFilter instanceof Rectangle || shapeFilter instanceof Point) || (nextFlexCell.cellStack.shapeFilterBoundingBox.getCrossesDateLine()) || rel == SpatialRelation.WITHIN)
          rel = null;
      }
      if (rel == null) {
        rel = nextCell.getShape().relate(shapeFilter);
      }
      return rel;
    }


    private void stopLevelIteration() {
      nextCellNumber = endCellNumber + 1;
    }

    //Populates into scratch the next cell in z-order TODO Hilbert ordering
    private boolean levelHasUntraversedCell() {
      if (nextCellNumber > endCellNumber) {
        nextCell = null;
        return false;
      }
      this.cell.cellStack.term.length = this.cell.cellLevel;
      //We must call this as we want the cell to invalidate its ShapeCache
      changeTailByte((byte) nextCellNumber);
      cell.reuse();
      ++nextCellNumber;
      return true;
    }
  }

  /**
   * A stack of flexCells with the following characteristics
   * - Lazy decoding of cells
   * - Cells from the same CellStack share BytesRef
   */
  private class CellStack {

    protected final FlexCell cells[];
    protected int lastDecodedLevel = 0;
    protected BytesRef term;

    //ShapeFilter bounding box and calculations
    private int shapeFilterXMin;
    private int shapeFilterXMax;
    private int shapeFilterYMin;
    private int shapeFilterYMax;
    private Rectangle shapeFilterBoundingBox;
    private Shape shapeFilter;

    public CellStack(int maxLevels, int xmin, int ymin) {
      this.cells = new FlexCell[maxLevels + 1];
      term = new BytesRef(maxLevels + 1); //+1 For leaf and this byteRef will be shared within the stack
      for (int level = maxLevels; level >= 0; --level) {
        if (level != maxLevels) {
          cells[level] = new FlexCell(cells[level + 1], this, level);
        } else {
          cells[level] = new FlexCell(null, this, level);
        }
      }
      //? The xmin,ymin needs to be set for the top cell. From there its decoded lazily a level at a time
      cells[0].setMinCornerCoordinates(xmin, ymin);
    }

    private void findIntegerBoundingBox(Shape shapeFilter) {
      if (shapeFilter != null && this.shapeFilter != shapeFilter) { // object equivalence?
        //TODO this remains same for a given FPT and given shape
        shapeFilterBoundingBox = shapeFilter.getBoundingBox();
        this.shapeFilterXMax = (int) ((shapeFilterBoundingBox.getMaxX() - bounds.getMinX()) * doubleToInt);
        this.shapeFilterXMin = (int) ((shapeFilterBoundingBox.getMinX() - bounds.getMinX()) * doubleToInt);
        this.shapeFilterYMax = (int) ((shapeFilterBoundingBox.getMaxY() - bounds.getMinY()) * doubleToInt);
        this.shapeFilterYMin = (int) ((shapeFilterBoundingBox.getMinY() - bounds.getMinY()) * doubleToInt);
        this.shapeFilter = shapeFilter;
      }

    }

    protected void decode(int cellLevel) {
      int xmin;
      int ymin;
      int row;
      int col;
      int c;
      int division;
      //decode all cells from the last decoded cell to the desired cell
      for (int i = lastDecodedLevel; i < cellLevel; i++) {
        xmin = cells[i].xMin;
        ymin = cells[i].yMin;
        c = term.bytes[term.offset + i] - 2;
        division = numberOfSubCellsAsExponentOfFour[i];
        col = (c >> division);
        row = (c - (1 << division) * col); // Is this worthwhile?
        xmin += gridSizes[i + 1] * col;
        ymin += gridSizes[i + 1] * row;
        cells[i + 1].setMinCornerCoordinates(xmin, ymin);
      }
      if (lastDecodedLevel < cellLevel) {
        lastDecodedLevel = cellLevel;
      }
    }


    /**
     * Invalidates the decoding of a cell forcing decoding to happen again
     *
     * @param cellLevel the Cell whose decoding is to be done
     */
    protected void invalidate(int cellLevel) {
      lastDecodedLevel = Math.max(cellLevel - 1, 0); //Note: Cell at level 0 is always decoded.
    }

  }
}
