import java.nio.FloatBuffer;

import java.util.*;

public class AxisAlignedTriangulation {
	
	public static class IntPosition2D {
		public int x, y;
		public IntPosition2D(int x, int y) {
			this.x = x;
			this.y = y;
		}

		public long lengthSquared() {
			return (long)x*x + (long)y*y;
		}
		
		public static long distSqr(IntPosition2D a, IntPosition2D b) {
			long dx = a.x - b.x;
			long dy = a.y - b.y;
			return dx*dx + dy*dy;
		}
		
		public static long dot(IntPosition2D a, IntPosition2D b) {
			return (long)a.x * b.x + (long)a.y * b.y;
		}
		
		public static IntPosition2D add(IntPosition2D a, IntPosition2D b, IntPosition2D dest) {
			if (dest == null)
				return new IntPosition2D(a.x + b.x, a.y + b.y);
			dest.set(a.x + b.x, a.y + b.y);
			return dest;
		}
		
		public static IntPosition2D sub(IntPosition2D a, IntPosition2D b, IntPosition2D dest) {
			if (dest == null)
				return new IntPosition2D(a.x - b.x, a.y - b.y);
			dest.set(a.x - b.x, a.y - b.y);
			return dest;
		}
		
		/** returns the signed triangle area * 2 of a triangle made of the points a, b, c in ccw order. */
		public static long triangleArea2(IntPosition2D a, IntPosition2D b, IntPosition2D c) {
			return (b.y - a.y) * (long)(c.x - a.x) - (b.x - a.x) * (long)(c.y - a.y);
		}
		
		/** "to" can be null */
		public void set(IntPosition2D to) {
			if (to != null) {
				x = to.x;
				y = to.y;
			}
			else {
				x = 0;
				y = 0;
			}
		}
		public void set(int x, int y) {
			this.x = x;
			this.y = y;
		}
		@Override
		public boolean equals(Object o) {
			IntPosition2D pos = (IntPosition2D)o;
			return x == pos.x && y == pos.y;
		}
		@Override
		public int hashCode() {
			return 667162477 * y + 119340589 * x;
		}
		@Override
		public String toString() {
			return String.format("%d, %d", x, y);
		}
	}
	
	private static void convexify(List<IntPosition2D> polyline, FloatBuffer positions, int start, int dir) {
		for (int i = start; i < polyline.size() - 1 && i > 0; i += dir) {
			IntPosition2D a = polyline.get(i - 1);
			IntPosition2D b = polyline.get(i);
			IntPosition2D c = polyline.get(i + 1);
			long area = IntPosition2D.triangleArea2(a, b, c);
			if (area >= 0) {
				polyline.remove(i);
				if (dir > 0)
					i--;
			}
			else
				break;
			if (area > 0) {
				positions.put(a.x);
				positions.put(a.y);
				positions.put(b.x);
				positions.put(b.y);
				positions.put(c.x);
				positions.put(c.y);
			}
		}
	}
	
	public static FloatBuffer triangulate(List<List<IntPosition2D>> polygons) {
		int triangles = 0;
		// remove all points on straight edges and count triangles
		List<List<IntPosition2D>> polygons_reduced = new ArrayList<>();
		for (List<IntPosition2D> polygon : polygons) {
			IntPosition2D last_last = polygon.get(polygon.size() - 2);
			IntPosition2D last = polygon.get(polygon.size() - 1);
			List<IntPosition2D> polygon_reduced = new ArrayList<>();
			for (int i = 0; i < polygon.size(); i++) {
				IntPosition2D cur = polygon.get(i);
				if (IntPosition2D.dot(IntPosition2D.sub(last_last, last, null), IntPosition2D.sub(last, cur, null)) == 0) {
					polygon_reduced.add(last);
				}
				last_last = last;
				last = cur;
			}
			polygons_reduced.add(polygon_reduced);
			triangles += polygon_reduced.size() - 2;
		}
		FloatBuffer positions = FloatBuffer.allocate(triangles * 3 * 2);
		
		for (List<IntPosition2D> polygon : polygons_reduced) {
			// triangulate using sweep line
			// here is the ways the sweep line can change in one step:
			// 1. it can gain a segment which is unconnected to anything
			// 2. it can loose a segment, which is inside of another segment convex polyline
			//   -> create all triangles to the polyline that have positive winding.
			// 3. a segment can be extended
			//   -> create triangles until the segment polyline is convex again.
			// 4. a segment can connect two neighboring segment convex polylines
			//   -> connect them up and add triangles until it is a convex polyline like in case 3
			
			assert polygon.size() % 2 == 0; // this is always true for axis aligned meshes without straight points
			Integer[] event_order = new Integer[polygon.size()];
			for (int i = 0; i < event_order.length; i++)
				event_order[i] = i;
			Arrays.parallelSort(event_order, (a, b) -> {
				IntPosition2D pa = polygon.get(a);
				IntPosition2D pb = polygon.get(b);
				if (pa.x != pb.x)
					return Integer.compare(pa.x, pb.x);
				// order the elements inside one sweepline event (same x) as well
				if (pa.y != pb.y)
					return Integer.compare(pa.y, pb.y);
				// this is the + crossing special case -> order by connection direction
				IntPosition2D pap = polygon.get((a + 1) % polygon.size());
				IntPosition2D pbp = polygon.get((b + 1) % polygon.size());
				IntPosition2D pan = polygon.get((a - 1 + polygon.size()) % polygon.size());
				IntPosition2D pbn = polygon.get((b - 1 + polygon.size()) % polygon.size());
				return Integer.compare(pap.y + pan.y, pbp.y + pbn.y);
			});
			
			// ordered list of segments (convex polylines) making up the sweep line
			List<List<IntPosition2D>> sweep_line = new ArrayList<>();
			for (int i = 0; i < event_order.length; i += 2) {
				// process one segment at a time
				int start_i = event_order[i];
				int end_i = event_order[i + 1];
				// these events are always (index-) neighbors in the polygon
				assert Math.abs(((start_i - end_i + 3*polygon.size()/2) % polygon.size()) - polygon.size()/2) == 1;
				IntPosition2D start = polygon.get(start_i);
				IntPosition2D end = polygon.get(end_i);
				assert start.x == end.x;
				assert start.y < end.y; // <= ensured by sorting, != ensured by data source
				// check if this edge gets added or removed from the sweep line
				boolean additive = (start_i - end_i + 1) % polygon.size() == 0;
				if (sweep_line.size() == 0) {
					assert additive;
					// add new segment to the split line
					List<IntPosition2D> new_segement = new ArrayList<>();
					new_segement.add(start);
					new_segement.add(end);
					sweep_line.add(new_segement);
					continue;
				}
				// find the segment in the sweep_line (HACK: use null as key and assume it is always passed in the same argument of the comparator)
				int sweep_line_i = Collections.binarySearch(sweep_line, null, (a, b) -> {
					assert b == null;
					return Integer.compare(a.get(a.size() - 1).y, start.y);
				});
				sweep_line_i = sweep_line_i >= 0 ? sweep_line_i : -sweep_line_i - 1;
				// check if the new segment is inside/connected to the old one in this place (this is the part up (with y <= start.y) to the new segment)
				List<IntPosition2D> convex_polyline = sweep_line_i < sweep_line.size() ? sweep_line.get(sweep_line_i) : null;
				if (convex_polyline != null && convex_polyline.get(0).y > end.y)
					convex_polyline = null;
				
				if (additive) {
					if (convex_polyline == null) {
						// add new segment to the split line
						List<IntPosition2D> new_segement = new ArrayList<>();
						new_segement.add(start);
						new_segement.add(end);
						sweep_line.add(sweep_line_i, new_segement);
					}
					else {
						// extend a segment to the negative side of this event
						boolean connect_negative = start.y == convex_polyline.get(convex_polyline.size() - 1).y;
						boolean connect_positive = end.y == convex_polyline.get(0).y;
						assert connect_negative || connect_positive;
						if (connect_negative) {
							convex_polyline.add(start);
							convex_polyline.add(end);
							convexify(convex_polyline, positions, convex_polyline.size() - 3, -1);
						}
						else { assert connect_positive;
							convex_polyline.add(0, end);
							convex_polyline.add(0, start);
							convexify(convex_polyline, positions, 2, 1);
						}
						
						// check if it's connected on the positive side of this event as well
						if (sweep_line_i + 1 < sweep_line.size() && sweep_line.get(sweep_line_i + 1).get(0).y == end.y) {
							assert connect_negative && !connect_positive;
							int index = convex_polyline.size();
							convex_polyline.addAll(sweep_line.get(sweep_line_i + 1));
							sweep_line.remove(sweep_line_i + 1);
							
							convexify(convex_polyline, positions, index, 1);
						}
					}
				}
				else {
					assert convex_polyline != null;
					// create new triangles from the connection to the convex part
					// -> convexify (remove points and return the position of the inserted segment)
					assert convex_polyline.get(0).y <= start.y && convex_polyline.get(convex_polyline.size() - 1).y >= end.y;
					// find an insert location which is visible from the point "start" and "end"
					// -> the most advanced vertex in the x direction is always visible from anything to its right
					int insert = 0;
					int best_x = Integer.MIN_VALUE;
					for (int j = 0; j < convex_polyline.size(); j++) {
						if (convex_polyline.get(j).x > best_x) {
							best_x = convex_polyline.get(j).x;
							insert = j;
						}
					}

					IntPosition2D insert_pos = convex_polyline.get(insert);
					// check area and only add triangles with positive area (yes even here)
					long area = IntPosition2D.triangleArea2(start, insert_pos, end);
					if (area > 0) {
						positions.put(start.x);
						positions.put(start.y);
						positions.put(insert_pos.x);
						positions.put(insert_pos.y);
						positions.put(end.x);
						positions.put(end.y);
					}

					List<IntPosition2D> convex_polyline1 = new ArrayList<>();
					convex_polyline1.addAll(convex_polyline.subList(0, insert + 1));
					convex_polyline1.add(start);
					
					List<IntPosition2D> convex_polyline2 = new ArrayList<>();
					convex_polyline2.add(end);
					convex_polyline2.addAll(convex_polyline.subList(insert, convex_polyline.size()));
					
					convexify(convex_polyline1, positions, convex_polyline1.size() - 2, -1);
					convexify(convex_polyline2, positions, 1, 1);
					
					sweep_line.remove(sweep_line_i);
					if (convex_polyline2.get(0).y < convex_polyline2.get(convex_polyline2.size() - 1).y)
						sweep_line.add(sweep_line_i, convex_polyline2);
					if (convex_polyline1.get(0).y < convex_polyline1.get(convex_polyline1.size() - 1).y)
						sweep_line.add(sweep_line_i, convex_polyline1);
				}
			}
		}
		positions.flip();
		return positions;
	}
}
