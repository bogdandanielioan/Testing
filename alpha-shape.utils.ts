import type { Outline, Position } from '../../../repair/models';
import { OutlineType } from '../../enums/outline-type.enum';
import { Shape, Vector2, Path } from 'three';
import { Delaunay } from 'd3-delaunay';

// Constants
const DEFAULT_SNAP_THRESHOLD = 0.2;
const DEFAULT_MIN_SEGMENT_LENGTH = 0.01;
const DEFAULT_ARC_SEGMENTS = 20;
const DEFAULT_ALPHA = 0.6;
const EPSILON = 1e-4;

/**
 * Creates alpha shapes (with holes) from all points in the provided outlines.
 */
export function createAlphaShapeFromOutlines(
  outlines: Outline[],
  alpha: number = DEFAULT_ALPHA
): Shape[] {
  // Collect and preprocess points
  const points = collectPointsFromOutlines(outlines);
  if (points.length < 3) return [];

  // Compute alpha shape boundaries
  const boundaries = computeAlphaShapeBoundaries(points, alpha);
  if (boundaries.length === 0) return [];

  // Convert to Three.js Shapes with proper holes
  return createShapesWithHoles(boundaries);
}

/**
 * Collects and deduplicates points from all outlines
 */
function collectPointsFromOutlines(outlines: Outline[]): Position[] {
  const allPoints: Position[] = [];
  outlines.forEach(outline => {
    if (outline.points) allPoints.push(...outline.points);
  });
  return deduplicatePoints(allPoints);
}

/**
 * Removes duplicate points (within epsilon tolerance)
 */
function deduplicatePoints(points: Position[], epsilon: number = EPSILON): Position[] {
  const seen = new Set<string>();
  return points.filter(p => {
    const key = `${Math.round(p.x / epsilon)},${Math.round(p.y / epsilon)}`;
    return seen.has(key) ? false : (seen.add(key), true);
  });
}

/**
 * Computes all boundaries (outer and holes) from the alpha shape
 */
function computeAlphaShapeBoundaries(points: Position[], alpha: number): Position[][] {
  if (points.length < 4) return [points];

  // Create Delaunay triangulation
  const delaunay = Delaunay.from(points.map(p => [p.x, p.y]));
  const triangles = extractTriangles(delaunay, points);

  // Find boundary edges based on alpha criterion
  const boundaryEdges = findAlphaBoundaryEdges(triangles, alpha);
  if (boundaryEdges.size === 0) return [];

  // Extract all closed boundaries from the edges
  const boundaries = extractClosedBoundaries(boundaryEdges);

  // Sort by area and determine outer boundary vs holes
  return sortAndClassifyBoundaries(boundaries);
}

/**
 * Extracts triangles from Delaunay triangulation
 */
function extractTriangles(delaunay: Delaunay<number>, points: Position[]): Triangle[] {
  const triangles: Triangle[] = [];
  for (let i = 0; i < delaunay.triangles.length; i += 3) {
    triangles.push({
      a: points[delaunay.triangles[i]],
      b: points[delaunay.triangles[i + 1]],
      c: points[delaunay.triangles[i + 2]]
    });
  }
  return triangles;
}

/**
 * Finds edges that form the alpha shape boundary
 */
function findAlphaBoundaryEdges(triangles: Triangle[], alpha: number): Map<string, Edge> {
  const edges = new Map<string, Edge>();

  for (const {a, b, c} of triangles) {
    const area = triangleArea(a, b, c);
    const maxEdge = Math.max(distance(a, b), distance(b, c), distance(c, a));
    const alphaMetric = area / (maxEdge ** 2);

    if (alphaMetric < alpha) {
      toggleEdge(a, b, edges);
      toggleEdge(b, c, edges);
      toggleEdge(c, a, edges);
    }
  }

  return edges;
}

/**
 * Toggles edge in the map (add if not present, remove if present)
 */
function toggleEdge(p1: Position, p2: Position, edges: Map<string, Edge>): void {
  const [a, b] = sortPoints(p1, p2);
  const key = `${a.x},${a.y}|${b.x},${b.y}`;
  edges.has(key) ? edges.delete(key) : edges.set(key, [p1, p2]);
}

/**
 * Extracts all closed boundaries from the edge map
 */
function extractClosedBoundaries(edges: Map<string, Edge>): Position[][] {
  const edgeList = Array.from(edges.values());
  const boundaries: Position[][] = [];

  while (edgeList.length > 0) {
    const boundary: Position[] = [edgeList[0][0], edgeList[0][1]];
    edgeList.splice(0, 1);

    let changed;
    do {
      changed = false;
      const last = boundary[boundary.length - 1];
      
      for (let i = 0; i < edgeList.length; i++) {
        const [a, b] = edgeList[i];
        if (pointsEqual(a, last)) {
          boundary.push(b);
          edgeList.splice(i, 1);
          changed = true;
          break;
        } else if (pointsEqual(b, last)) {
          boundary.push(a);
          edgeList.splice(i, 1);
          changed = true;
          break;
        }
      }
    } while (changed && edgeList.length > 0);

    // Close the loop if start and end points match
    if (pointsEqual(boundary[0], boundary[boundary.length - 1])) {
      boundaries.push(boundary);
    }
  }

  return boundaries;
}

/**
 * Sorts boundaries by area and classifies them as outer or holes
 */
function sortAndClassifyBoundaries(boundaries: Position[][]): Position[][] {
  if (boundaries.length === 0) return [];
  
  // Calculate signed area for each boundary
  const withAreas = boundaries.map(b => ({
    boundary: b,
    area: polygonSignedArea(b)
  }));

  // Sort by absolute area (descending)
  withAreas.sort((a, b) => Math.abs(b.area) - Math.abs(a.area));

  // The largest boundary is the outer shape
  const result = [withAreas[0].boundary];

  // Other boundaries are holes if they're inside the outer boundary
  const outer = withAreas[0].boundary;
  for (let i = 1; i < withAreas.length; i++) {
    const boundary = withAreas[i].boundary;
    if (isPointInPolygon(boundary[0], outer)) {
      // Reverse hole boundaries to have opposite winding
      result.push(withAreas[i].area > 0 ? boundary.reverse() : boundary);
    }
  }

  return result;
}

/**
 * Creates Three.js Shapes with holes from boundaries
 */
function createShapesWithHoles(boundaries: Position[][]): Shape[] {
  if (boundaries.length === 0) return [];

  const shapes: Shape[] = [];
  const outerBoundary = boundaries[0];
  const shape = new Shape(outerBoundary.map(p => new Vector2(p.x, p.y)));

  // Add holes (remaining boundaries)
  for (let i = 1; i < boundaries.length; i++) {
    const hole = boundaries[i];
    shape.holes.push(new Path(hole.map(p => new Vector2(p.x, p.y))));
  }

  shapes.push(shape);
  return shapes;
}

/**
 * Line segment preprocessing functions
 */
export function preprocessLineSegments(
  lines: Outline[],
  snapThreshold: number = DEFAULT_SNAP_THRESHOLD
): Outline[] {
  // Filter valid lines
  let result = lines.filter(o => 
    o.type === OutlineType.Line && 
    o.points?.length === 2 &&
    distance(o.points[0], o.points[1]) > DEFAULT_MIN_SEGMENT_LENGTH
  );

  // Process lines
  result = removeDuplicateLines(result, snapThreshold);
  snapLineEndpoints(result, snapThreshold);
  connectCloseEnds(result, snapThreshold);

  return result;
}

function removeDuplicateLines(lines: Outline[], threshold: number): Outline[] {
  const seen = new Set<string>();
  return lines.filter(line => {
    const [a, b] = sortPoints(line.points[0], line.points[1]);
    const key = `${roundValue(a.x, threshold)}_${roundValue(a.y, threshold)}|${
      roundValue(b.x, threshold)}_${roundValue(b.y, threshold)}`;
    return seen.has(key) ? false : (seen.add(key), true);
  });
}

function snapLineEndpoints(lines: Outline[], threshold: number): void {
  const allPoints = lines.flatMap(line => line.points);
  
  for (let i = 0; i < allPoints.length; i++) {
    for (let j = i + 1; j < allPoints.length; j++) {
      if (arePointsClose(allPoints[i], allPoints[j], threshold)) {
        allPoints[j].x = allPoints[i].x;
        allPoints[j].y = allPoints[i].y;
      }
    }
  }
}

function connectCloseEnds(lines: Outline[], threshold: number): void {
  for (let i = 0; i < lines.length; i++) {
    const endA = lines[i].points[1];
    for (let j = 0; j < lines.length; j++) {
      if (i === j) continue;
      const startB = lines[j].points[0];
      if (arePointsClose(endA, startB, threshold)) {
        startB.x = endA.x;
        startB.y = endA.y;
      }
    }
  }
}

/**
 * Arc and circle conversion functions
 */
export function arcToLineSegments(
  arc: Outline, 
  segments: number = DEFAULT_ARC_SEGMENTS
): Outline | null {
  if (!arc.start || !arc.end || arc.startTheta == null || arc.endTheta == null) return null;
  
  const points: Position[] = [];
  const startAngle = (arc.startTheta / 10) * (Math.PI / 180);
  const endAngle = (arc.endTheta / 10) * (Math.PI / 180);
  const angleDelta = endAngle - startAngle;
  const center = {
    x: (arc.start.x + arc.end.x) / 2,
    y: (arc.start.y + arc.end.y) / 2
  };
  const radiusX = Math.abs(arc.start.x - arc.end.x) / 2;
  const radiusY = Math.abs(arc.start.y - arc.end.y) / 2;

  for (let i = 0; i <= segments; i++) {
    const t = startAngle + angleDelta * (i / segments);
    points.push({
      x: center.x + radiusX * Math.cos(t),
      y: center.y + radiusY * Math.sin(t)
    });
  }

  return { type: OutlineType.Polygon, points };
}

export function circleToLineSegments(
  circle: Outline,
  segments: number = DEFAULT_ARC_SEGMENTS
): Outline | null {
  if (!circle.location || circle.radius == null) return null;
  
  const points: Position[] = [];
  for (let i = 0; i <= segments; i++) {
    const angle = (2 * Math.PI * i) / segments;
    points.push({
      x: circle.location.x + circle.radius * Math.cos(angle),
      y: circle.location.y + circle.radius * Math.sin(angle)
    });
  }

  return { type: OutlineType.Polygon, points };
}

/**
 * Geometric utility functions
 */
function distance(a: Position, b: Position): number {
  const dx = a.x - b.x;
  const dy = a.y - b.y;
  return Math.sqrt(dx * dx + dy * dy);
}

function pointsEqual(a: Position, b: Position, epsilon: number = EPSILON): boolean {
  return Math.abs(a.x - b.x) < epsilon && Math.abs(a.y - b.y) < epsilon;
}

function sortPoints(a: Position, b: Position): [Position, Position] {
  return a.x < b.x || (a.x === b.x && a.y < b.y) ? [a, b] : [b, a];
}

function arePointsClose(a: Position, b: Position, threshold: number): boolean {
  return distance(a, b) <= threshold;
}

function roundValue(val: number, inc: number): number {
  return Math.round(val / inc) * inc;
}

function triangleArea(a: Position, b: Position, c: Position): number {
  return Math.abs(
    (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y)) / 2
  );
}

function polygonSignedArea(points: Position[]): number {
  let area = 0;
  const n = points.length;
  for (let i = 0, j = n - 1; i < n; j = i++) {
    area += (points[j].x + points[i].x) * (points[j].y - points[i].y);
  }
  return area / 2;
}

function isPointInPolygon(point: Position, polygon: Position[]): boolean {
  let inside = false;
  const x = point.x, y = point.y;
  for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
    const xi = polygon[i].x, yi = polygon[i].y;
    const xj = polygon[j].x, yj = polygon[j].y;

    const intersect = ((yi > y) !== (yj > y)) &&
      (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
    if (intersect) inside = !inside;
  }
  return inside;
}

// Type definitions
type Edge = [Position, Position];
interface Triangle {
  a: Position;
  b: Position;
  c: Position;
}
