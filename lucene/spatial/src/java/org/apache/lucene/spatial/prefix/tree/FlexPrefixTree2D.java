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

import java.io.PrintStream;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Locale;

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
  private final double xmin;
  private final double xmax;
  private final double ymin;
  private final double ymax;
  private final double xmid;
  private final double ymid;

  private final double gridW;
  public final double gridH;

  final double[] levelW;
  final double[] levelH;
  final int[]    levelS; // side
  final int[]    levelN; // number

  public FlexPrefixTree2D(SpatialContext ctx, Rectangle bounds, int maxLevels) {
    super(ctx, maxLevels);
    this.xmin = bounds.getMinX();
    this.xmax = bounds.getMaxX();
    this.ymin = bounds.getMinY();
    this.ymax = bounds.getMaxY();

    levelW = new double[maxLevels];
    levelH = new double[maxLevels];
    levelS = new int[maxLevels];
    levelN = new int[maxLevels];

    gridW = xmax - xmin;
    gridH = ymax - ymin;
    this.xmid = xmin + gridW/2.0;
    this.ymid = ymin + gridH/2.0;
    levelW[0] = gridW/2.0;
    levelH[0] = gridH/2.0;
    levelS[0] = 2;
    levelN[0] = 4;

    for (int i = 1; i < levelW.length; i++) {
      levelW[i] = levelW[i - 1] / 2.0;
      levelH[i] = levelH[i - 1] / 2.0;
      levelS[i] = levelS[i - 1] * 2;
      levelN[i] = levelN[i - 1] * 4;
    }
  }

  public FlexPrefixTree2D(SpatialContext ctx) {
    this(ctx, DEFAULT_MAX_LEVELS);
  }

  public FlexPrefixTree2D(SpatialContext ctx, int maxLevels) {
    this(ctx, ctx.getWorldBounds(), maxLevels);
  }

  @Override
  public Cell getWorldCell() {
    return newCellStack(maxLevels)[0];
  }

  protected FlexCell[] newCellStack(int levels) {
    levels++;
    final FlexCell[] cellsByLevel = new FlexCell[levels + 1];
    final BytesRef term = new BytesRef(levels+1); //For leaf
    for (int level = 0; level <= levels; level++) {
      cellsByLevel[level] = new FlexCell(cellsByLevel,term.bytes,term.offset,level);
    }
    return cellsByLevel;
  }

  public void printInfo(PrintStream out) {
    NumberFormat nf = NumberFormat.getNumberInstance(Locale.ROOT);
    nf.setMaximumFractionDigits(5);
    nf.setMinimumFractionDigits(5);
    nf.setMinimumIntegerDigits(3);

    for (int i = 0; i < maxLevels; i++) {
      out.println(i + "]\t" + nf.format(levelW[i]) + "\t" + nf.format(levelH[i]) + "\t" +
          levelS[i] + "\t" + (levelS[i] * levelS[i]));
    }
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
    List<Cell> cells = new ArrayList<>(1);
    build(xmid, ymid, 0, cells, new BytesRef(maxLevels+1), ctx.makePoint(p.getX(),p.getY()), level,cellsByLevel);
    return cells.get(0);//note cells could be longer if p on edge
  }

  private void build(double x,double y,int level,List<Cell> matches,BytesRef str,Shape shape,int maxLevel,FlexCell[] cellsByLevel) {
    assert str.length == level;
    double w = levelW[level] / 2;
    double h = levelH[level] / 2;

    // Z-Order
    // http://en.wikipedia.org/wiki/Z-order_%28curve%29
    checkBattenberg((byte)0x02, x - w, y + h, level, matches, str, shape, maxLevel,cellsByLevel);
    checkBattenberg((byte)0x03, x + w, y + h, level, matches, str, shape, maxLevel,cellsByLevel);
    checkBattenberg((byte)0x04, x - w, y - h, level, matches, str, shape, maxLevel,cellsByLevel);
    checkBattenberg((byte)0x05, x + w, y - h, level, matches, str, shape, maxLevel,cellsByLevel);

    // possibly consider hilbert curve
    // http://en.wikipedia.org/wiki/Hilbert_curve
    // http://blog.notdot.net/2009/11/Damn-Cool-Algorithms-Spatial-indexing-with-Quadtrees-and-Hilbert-Curves
    // if we actually use the range property in the query, this could be useful
  }

  private void checkBattenberg(byte c,double cx,double cy,int level,List<Cell> matches,BytesRef str,Shape shape,int maxLevel,FlexCell[] cellsByLevel) {
    assert str.length == level;
    assert str.offset == 0;
    double w = levelW[level] / 2;
    double h = levelH[level] / 2;

    int strlen = str.length;
    Rectangle rectangle = ctx.makeRectangle(cx - w, cx + w, cy - h, cy + h);
    SpatialRelation v = shape.relate(rectangle);
    if (SpatialRelation.CONTAINS == v) {
      str.bytes[str.length++] = c;//append
      //str.append(SpatialPrefixGrid.COVER);
      matches.add(cellsByLevel[level].reuse(BytesRef.deepCopyOf(str), v.transpose()));
    } else if (SpatialRelation.DISJOINT == v) {
      // nothing
    } else { // SpatialRelation.WITHIN, SpatialRelation.INTERSECTS
      str.bytes[str.length++] = c;//append

      int nextLevel = level+1;
      if (nextLevel >= maxLevel) {
        //str.append(SpatialPrefixGrid.INTERSECTS);
        matches.add(cellsByLevel[level].reuse(BytesRef.deepCopyOf(str), v.transpose()));
      } else {
        build(cx, cy, nextLevel, matches, str, shape, maxLevel,cellsByLevel);
      }
    }
    str.length = strlen;
  }

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

  @Override
  public CellIterator getTreeCellIterator(Shape shape, int detailLevel) {
    if (!(shape instanceof Point))
      return super.getTreeCellIterator(shape, detailLevel);

    //This specialization is here because the legacy implementations don't have a fast implementation of
    // cell.getSubCells(point). It's fastest here to encode the full bytes for detailLevel, and create
    // subcells from the bytesRef in a loop. This avoids an O(N^2) encode, and we have O(N) instead.

    FlexCell worldCell = (FlexCell)getWorldCell();
    Cell cell = getCell((Point) shape, detailLevel,worldCell.cellsByLevel);
    assert !cell.isLeaf() && cell instanceof FlexCell;
    BytesRef fullBytes = cell.getTokenBytesNoLeaf(null);
    //fill in reverse order to be sorted
    Cell[] cells = new Cell[detailLevel];
    for (int i = 1; i < detailLevel; i++) {
      fullBytes.length = i;
      Cell parentCell = readCell(fullBytes, null);
      cells[i-1] = parentCell;
    }
    cells[detailLevel-1] = cell;
    return new FilterCellIterator(Arrays.asList(cells).iterator(), null);//null filter
  }


  public class FlexCell implements Cell {


    //Arguably we could simply use a BytesRef, using an extra Object.
    private byte[] bytes;//generally bigger to potentially hold a leaf
    private int b_off;
    private int b_len;//doesn't reflect leaf; same as getLevel()

    protected boolean isLeaf;

    /**
     * When set via getSubCells(filter), it is the relationship between this cell
     * and the given shape filter. Doesn't participate in shape equality.
     */
    protected SpatialRelation shapeRel;

    protected Shape shape;//cached

    final FlexCell[] cellsByLevel;

    final FlexPrefixTreeIterator cellIterator;

    /** Warning: Refers to the same bytes (no copy). If {@link #setLeaf()} is subsequently called then it
     * may modify bytes. */
    FlexCell(FlexCell[] cellsByLevel,byte[] bytes, int off, int len) {
      super();
      this.cellIterator = new FlexPrefixTreeIterator();
      this.cellsByLevel = cellsByLevel;
      this.bytes = bytes;
      this.b_off = off;
      this.b_len = len;
      readLeafAdjust();
    }

    FlexCell(FlexCell[] cellsByLevel,BytesRef str, SpatialRelation shapeRel) {
      this(cellsByLevel,str.bytes, str.offset, str.length);
      this.shapeRel = shapeRel;
    }

    private Cell reuse(byte[] bytes, int off, int len) {
      this.shape =null;
      this.isLeaf = false;
      this.shapeRel = null;
      this.bytes = bytes;
      this.b_off = off;
      this.b_len = len;
      readLeafAdjust();
      return this;
    }


    protected Cell reuse(BytesRef str, SpatialRelation shapeRel){
      this.reuse(str.bytes, str.offset, str.length);
      this.shapeRel = shapeRel;
      return this;
    }

    protected FlexPrefixTree2D getGrid() { return FlexPrefixTree2D.this; }

    /**
     * Gets the cells at the next grid cell level that covers this cell.
     * Precondition: Never called when getLevel() == maxLevel.
     *
     * @return A set of cells (no dups), sorted, modifiable, not empty, not null.
     */

    public int getSubCellsSize() {
      return 4;
    }


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
      if (shape == null)
        shape = makeShape();
      return shape;
    }

    private Rectangle makeShape() {
      BytesRef token = getTokenBytesNoLeaf(null);
      double xmin = FlexPrefixTree2D.this.xmin;
      double ymin = FlexPrefixTree2D.this.ymin;

      for (int i = 0; i < token.length; i++) {
        byte c = token.bytes[token.offset + i];
        switch (c) {
          case 0x02:
            ymin += levelH[i];
            break;
          case 0x03:
            xmin += levelW[i];
            ymin += levelH[i];
            break;
          case 0x04:
            break;//nothing really
          case 0x05:
            xmin += levelW[i];
            break;
          default:
            throw new RuntimeException("unexpected char: " + c);
        }
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
      return ctx.makeRectangle(xmin, xmin + width, ymin, ymin + height);
    }

    private static final byte LEAF_BYTE = 0x01;//NOTE: must sort before letters & numbers

    protected void readCell(BytesRef bytes) {
      shapeRel = null;
      shape = null;
      this.bytes = bytes.bytes;
      this.b_off = bytes.offset;
      this.b_len = bytes.length;
      readLeafAdjust();
    }

    private void readLeafAdjust() {
      isLeaf = (b_len > 0 && bytes[b_off + b_len - 1] == LEAF_BYTE);
      if (isLeaf)
        b_len--;
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
      if (result.bytes.length < result.offset + result.length + 1) {
        assert false : "Not supposed to happen; performance bug";
        byte[] copy = new byte[result.length + 1];
        System.arraycopy(result.bytes, result.offset, copy, 0, result.length - 1);
        result.bytes = copy;
        result.offset = 0;
      }
      result.bytes[result.offset + result.length++] = LEAF_BYTE;
      return result;
    }

    @Override
    public BytesRef getTokenBytesNoLeaf(BytesRef result) {
      if (result == null)
        return new BytesRef(bytes, b_off, b_len);
      result.bytes = bytes;
      result.offset = b_off;
      result.length = b_len;
      return result;
    }

    @Override
    public int getLevel() {
      return b_len;
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
        BytesRef target = new BytesRef();
        return cellIterator.initIter(source, target, cellsByLevel[getLevel() + 1], shapeFilter, 0x02, 0x05);
      }
    }

    @Override
    public boolean isPrefixOf(Cell c) {
      //Note: this only works when each level uses a whole number of bytes.
      FlexCell cell = (FlexCell)c;
      boolean result = sliceEquals(cell.bytes, cell.b_off, cell.b_len, bytes, b_off, b_len);
      assert result == StringHelper.startsWith(c.getTokenBytesNoLeaf(null), getTokenBytesNoLeaf(null));
      return result;
    }

    /** Copied from {@link org.apache.lucene.util.StringHelper#startsWith(BytesRef, BytesRef)}
     *  which calls this. This is to avoid creating a BytesRef.  */
    private boolean sliceEquals(byte[] sliceToTest_bytes, int sliceToTest_offset, int sliceToTest_length,
                                       byte[] other_bytes, int other_offset, int other_length) {
      if (sliceToTest_length < other_length) {
        return false;
      }
      int i = sliceToTest_offset;
      int j = other_offset;
      final int k = other_offset + other_length;

      while (j < k) {
        if (sliceToTest_bytes[i++] != other_bytes[j++]) {
          return false;
        }
      }

      return true;
    }

    @Override
    public int compareToNoLeaf(Cell fromCell) {
      FlexCell b = (FlexCell) fromCell;
      return compare(bytes, b_off, b_len, b.bytes, b.b_off, b.b_len);
    }

    /** Copied from {@link BytesRef#compareTo(BytesRef)}.
     * This is to avoid creating a BytesRef. */
    private int compare(byte[] aBytes, int aUpto, int a_length, byte[] bBytes, int bUpto, int b_length) {
      final int aStop = aUpto + Math.min(a_length, b_length);
      while(aUpto < aStop) {
        int aByte = aBytes[aUpto++] & 0xff;
        int bByte = bBytes[bUpto++] & 0xff;

        int diff = aByte - bByte;
        if (diff != 0) {
          return diff;
        }
      }

      // One is a prefix of the other, or, they are equal:
      return a_length - b_length;
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

    /*@Override
    public int hashCode() {
      return getTokenBytesWithLeaf(null).hashCode();
    }*/

    @Override
    public String toString() {
      //this method isn't "normally" called; just in asserts/tests
      return getTokenBytesWithLeaf(null).utf8ToString();
      }

  }//QuadCell

  /**
   * An Iterator for FlexCells. This iterator reuses cells at a level and iterates over the siblings
   */
  private class FlexPrefixTreeIterator extends CellIterator{

    private BytesRef source;
    private BytesRef target;
    private int start;
    private int end;
    private FlexCell scratch;
    private Shape shapeFilter;

    protected CellIterator initIter(BytesRef source,BytesRef target,FlexCell scratch,Shape shapeFilter,int start,int end){
      this.nextCell = null;
      this.thisCell = null;
      this.source = source;
      this.target = target;
      this.scratch = scratch;
      this.shapeFilter = shapeFilter;
      this.start = start;
      this.end = end;
      return this;
    }

    private BytesRef concat(BytesRef source, byte b, BytesRef target) {
      assert target.offset == 0;
      target.bytes = new byte[source.length + 2];//+2 for new char + potential leaf
      target.length = 0;
      target.copyBytes(source);
      target.bytes[target.length++] = b;
      return target;
    }

    private boolean hasNextCell(){
      //Iterate over the cells
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
              nextCell.setLeaf();
            return true;
          }
        }
      }
      return false;
    }
  }
}
