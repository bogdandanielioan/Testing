import type { Outline, Position } from '../../../repair/models';
import { OutlineType } from '../../enums/outline-type.enum';
import { Shape, Vector2 } from 'three';
import { Delaunay } from 'd3-delaunay';

const DEFAULT_SNAP_THRESHOLD = 0.2;
const DEFAULT_MIN_SEGMENT_LENGTH = 0.01;
const DEFAULT_ARC_SEGMENTS = 20;

/**
 * Creates an alpha shape (a single Shape) from all points in the provided outlines.
 */
export function createAlphaShapeFromOutlines(
  outlines: Outline[],
  alpha: number = 0.6
): Shape[] {
  // Collect & deduplicate all points
  const allPoints: Position[] = [];
  outlines.forEach(outline => outline.points && allPoints.push(...outline.points));
  const deduped = deduplicatePoints(allPoints);

  if (deduped.length < 3) return [];
  const boundary = computeAlphaShape(deduped, alpha);
  if (!boundary.length) return [];

  return [new Shape(boundary.map(p => new Vector2(p.x, p.y)))];
}

function deduplicatePoints(points: Position[], epsilon = 1e-4): Position[] {
  const seen = new Set<string>();

  return points.filter(p => {
    const key = `${Math.round(p.x / epsilon)},${Math.round(p.y / epsilon)}`;
    if (seen.has(key)) return false;
    seen.add(key);

    return true;
  });
}

function computeAlphaShape(points: Position[], alpha: number): Position[] {
  if (points.length < 4) return points;

  const coords = points.map(p => [p.x, p.y]) as [number, number][];
  const delaunay = Delaunay.from(coords);

  // Gather triangles
  const triangles: [number, number, number][] = [];
  for (let i = 0; i < delaunay.triangles.length; i += 3) {
    triangles.push([
      delaunay.triangles[i],
      delaunay.triangles[i + 1],
      delaunay.triangles[i + 2],
    ]);
  }

  // Compute edges that form the alpha shape boundary
  const edges = new Map<string, [Position, Position]>();
  for (const [i1, i2, i3] of triangles) {
    const a = points[i1], b = points[i2], c = points[i3];
    const area = triangleArea(a, b, c);
    const maxEdge = Math.max(
      distance(a, b), distance(b, c), distance(c, a)
    );
    const alphaMetric = area / (maxEdge ** 2);
    if (alphaMetric < alpha) {
      addOrRemoveEdge(a, b, edges);
      addOrRemoveEdge(b, c, edges);
      addOrRemoveEdge(c, a, edges);
    }
  }

  // Reconstruct boundary by ordering the edges
  const edgeList = Array.from(edges.values());
  if (!edgeList.length) return [];

  const boundary: Position[] = [edgeList[0][0], edgeList[0][1]];
  edgeList.splice(0, 1);

  while (edgeList.length) {
    const last = boundary[boundary.length - 1];
    const idx = edgeList.findIndex(([p1, p2]) =>
      (p1.x === last.x && p1.y === last.y) ||
        (p2.x === last.x && p2.y === last.y)
    );
    if (idx === -1) break;
    const [a, b] = edgeList.splice(idx, 1)[0];
    boundary.push((a.x === last.x && a.y === last.y) ? b : a);
  }

  return boundary;
}

function triangleArea(a: Position, b: Position, c: Position): number {
  const ab = distance(a, b);
  const bc = distance(b, c);
  const ca = distance(c, a);
  const s = (ab + bc + ca) / 2;

  return Math.sqrt(s * (s - ab) * (s - bc) * (s - ca));
}

function addOrRemoveEdge(p1: Position, p2: Position, edges: Map<string, [Position, Position]>) {
  // Sort endpoints for a canonical key
  const [ptA, ptB] = [p1, p2].sort((a, b) => (a.x !== b.x ? a.x - b.x : a.y - b.y));
  const key = `${ptA.x},${ptA.y}|${ptB.x},${ptB.y}`;
  edges.has(key) ? edges.delete(key) : edges.set(key, [p1, p2]);
}

/**
 * Preprocesses line segments: filters, deduplicates, snaps, etc.
 */
export function preprocessLineSegments(
  lines: Outline[],
  snapThreshold: number = DEFAULT_SNAP_THRESHOLD
): Outline[] {
  let result = lines.filter(o => o.type === OutlineType.Line && o.points?.length === 2);
  result = removeDuplicateLines(result, snapThreshold);
  result = result.filter(line => distance(line.points[0], line.points[1]) > DEFAULT_MIN_SEGMENT_LENGTH);
  snapLineEndpoints(result, snapThreshold);
  connectCloseEnds(result, snapThreshold);

  return result;
}

function removeDuplicateLines(lines: Outline[], snapThreshold: number): Outline[] {
  const seen = new Set<string>();
  const output: Outline[] = [];

  for (const line of lines) {
    const [p1, p2] = line.points;
    const [a, b] = sortPoints(p1, p2);
    const key = `${roundValue(a.x, snapThreshold)}_${roundValue(a.y, snapThreshold)}|` +
        `${roundValue(b.x, snapThreshold)}_${roundValue(b.y, snapThreshold)}`;

    if (!seen.has(key)) {
      seen.add(key);
      output.push(line);
    }
  }

  return output;
}

function snapLineEndpoints(lines: Outline[], threshold: number): void {
  const allPoints: Position[] = [];
  lines.forEach(line => {
    allPoints.push(line.points[0], line.points[1]);
  });

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
    const [startA, endA] = lines[i].points;
    for (let j = i + 1; j < lines.length; j++) {
      const [startB, endB] = lines[j].points;
      if (arePointsClose(endA, startB, threshold)) {
        startB.x = endA.x;
        startB.y = endA.y;
      }
    }
  }
}

function arcToLineSegments(arc: Outline, segments: number = DEFAULT_ARC_SEGMENTS): Outline | null {
  if (!arc.start || !arc.end || arc.startTheta == null || arc.endTheta == null) return null;
  const points: Position[] = [];

  const start = (arc.startTheta / 10) * (Math.PI / 180);
  const end   = (arc.endTheta   / 10) * (Math.PI / 180);
  const cx = (arc.start.x + arc.end.x) / 2;
  const cy = (arc.start.y + arc.end.y) / 2;
  const rx = Math.abs(arc.start.x - arc.end.x) / 2;
  const ry = Math.abs(arc.start.y - arc.end.y) / 2;
  const sweep = end - start;

  for (let i = 0; i <= segments; i++) {
    const t = start + sweep * (i / segments);
    points.push({ x: cx + rx * Math.cos(t), y: cy + ry * Math.sin(t) });
  }

  return { type: OutlineType.Polygon, points };
}

function circleToLineSegments(circle: Outline, segments: number = DEFAULT_ARC_SEGMENTS): Outline | null {
  if (!circle.location || circle.radius == null) return null;
  const points: Position[] = [];

  for (let i = 0; i <= segments; i++) {
    const angle = (2 * Math.PI * i) / segments;
    points.push({
      x: circle.location.x + circle.radius * Math.cos(angle),
      y: circle.location.y + circle.radius * Math.sin(angle),
    });
  }

  return { type: OutlineType.Polygon, points };
}

/* Utility functions */
function distance(a: Position, b: Position): number {
  const dx = a.x - b.x;
  const dy = a.y - b.y;

  return Math.sqrt(dx * dx + dy * dy);
}

function sortPoints(a: Position, b: Position): [Position, Position] {
  return (a.x < b.x || (a.x === b.x && a.y <= b.y)) ? [a, b] : [b, a];
}

function arePointsClose(a: Position, b: Position, threshold: number): boolean {
  const dx = a.x - b.x;
  const dy = a.y - b.y;

  return dx * dx + dy * dy <= threshold * threshold;
}

function roundValue(val: number, inc: number): number {
  return Math.round(val / inc) * inc;
}

// Export arc/circle conversion if needed
export { arcToLineSegments, circleToLineSegments };
