module msdf.shape;
import msdf.contour;
import msdf.scanline;
import msdf.segment;
import inmath;
import std.algorithm.sorting;

enum MSDFGEN_CORNER_DOT_EPSILON = 0.000001;
enum MSDFGEN_DECONVERGENCE_FACTOR = 0.000001;

private {
    void _msdfDeconvergeEdge(ref EdgeSegment edge, int param) {
        {
            if (QuadraticSegment quadraticSegment = cast(QuadraticSegment)edge)
                edge = quadraticSegment.convertToCubic();
        }
        {
            if (CubicSegment quadraticSegment = cast(CubicSegment)edge)
                quadraticSegment.deconverge(param, MSDFGEN_DECONVERGENCE_FACTOR);
        }
    }
}

/// Vector shape representation.
class Shape {
public:
    struct Bounds {
        double l, b, r, t;
    }

    /// The list of contours the shape consists of.
    Contour[] contours;
    /// Specifies whether the shape uses bottom-to-top (false) or top-to-bottom (true) Y coordinates.
    bool inverseYAxis;

    this() {
        this.inverseYAxis = false;
    }

    /// Adds a contour.
    void addContour(Contour contour) {
        this.contours ~= contour;
    }
    
    /// Normalizes the shape geometry for distance field generation.
    void normalize() {
        foreach(ref contour; contours) {
            if (contour.edges.length == 1) {
                EdgeSegment[3] parts;
                contour.edges[0].splitInThirds(parts[0], parts[1], parts[2]);
                contour.edges.length = 0;
                contour.edges ~= parts[0];
                contour.edges ~= parts[1];
                contour.edges ~= parts[2];
            } else {
                EdgeSegment prevEdge = contour.edges[$-1];
                foreach(ref edge; contour.edges) {
                    vec2d prevDir = prevEdge.direction(1).normalized;
                    vec2d curDir = edge.direction(0).normalized;
                    if (dot(prevDir, curDir) < MSDFGEN_CORNER_DOT_EPSILON) {
                        _msdfDeconvergeEdge(prevEdge, 1);
                        _msdfDeconvergeEdge(edge, 0);
                    }
                    prevEdge = edge;
                }
            }
        }
    }

    /// Performs basic checks to determine if the object represents a valid shape.
    bool validate() const {
        foreach(ref contour; contours) {
            if (contour.edges.length != 0) {
                vec2d corner = contour.edges[$-1].point(1);
                foreach(ref edge; contour.edges) {
                    if (!edge) return false;
                    if (edge.point(0) != corner) return false;
                    corner = edge.point(1);
                }
            }
        }
        return true;
    }

    /// Adjusts the bounding box to fit the shape.
    void bound(ref double l, ref double b, ref double r, ref double t) const {
        foreach(ref contour; contours) 
            contour.bound(l, b, r, t);
    }

    /// Adjusts the bounding box to fit the shape border's mitered corners.
    void boundMiters(ref double l, ref double b, ref double r, ref double t, double border, double miterLimit, int polarity) const {
        foreach(ref contour; contours) 
            contour.boundMiters(l, b, r, t, border, miterLimit, polarity);
    }

    /// Computes the minimum bounding box that fits the shape, optionally with a (mitered) border.
    Bounds getBounds(double border = 0, double miterLimit = 0, int polarity = 0) const {
        static const double LARGE_VALUE = 1e240;
        Bounds bounds = { +LARGE_VALUE, +LARGE_VALUE, -LARGE_VALUE, -LARGE_VALUE };
        bound(bounds.l, bounds.b, bounds.r, bounds.t);
        if (border > 0) {
            bounds.l -= border, bounds.b -= border;
            bounds.r += border, bounds.t += border;
            if (miterLimit > 0)
                boundMiters(bounds.l, bounds.b, bounds.r, bounds.t, border, miterLimit, polarity);
        }
        return bounds;
    }

    /// Outputs the scanline that intersects the shape at y.
    void scanline(ref Scanline line, double y) const {
        Scanline.Intersection[] intersections;
        double[3] x;
        int[3] dy;
        foreach (ref contour; contours) {
            foreach (ref edge; contour.edges) {
                int n = edge.scanlineIntersections(x, dy, y);
                for (int i = 0; i < n; ++i) {
                    Scanline.Intersection intersection = { x[i], dy[i] };
                    intersections ~= intersection;
                }
            }
        }
        line.setIntersections(intersections);
    }

    /// Returns the total number of edge segments
    int edgeCount() const {
        int total = 0;
        foreach(ref contour; contours) total += cast(int)contour.edges.length;
        return total;
    }

    /// Assumes its contours are unoriented (even-odd fill rule). Attempts to orient them to conform to the non-zero winding rule.
    void orientContours() {
        struct Intersection {
            double x;
            int direction;
            int contourIndex;
        }

        const double ratio = .5*(sqrt(5.0)-1.0); // an irrational number to minimize chance of intersecting a corner or other point of interest
        int[] orientations = new int[contours.length];
        Intersection[] intersections;
        for (int i = 0; i < cast(int)contours.length; ++i) {
            if (!orientations[i] && contours[i].edges.length != 0) {
                // Find an Y that crosses the contour
                double y0 = contours[i].edges[$-1].point(0).y;
                double y1 = y0;
                foreach (ref edge; contours[i].edges)
                    y1 = edge.point(1).y;
                foreach (ref edge; contours[i].edges)
                    y1 = edge.point(ratio).y; // in case all endpoints are in a horizontal line
                double y = mix(y0, y1, ratio);
                // Scanline through whole shape at Y
                double[3] x;
                int[3] dy;
                foreach(j; 0..contours.length) {
                    foreach (ref edge; contours[j].edges) {
                        int n = edge.scanlineIntersections(x, dy, y);
                        for (int k = 0; k < n; ++k) {
                            Intersection intersection = { x[k], dy[k], j };
                            intersections ~= intersection;
                        }
                    }
                }
                
                sort!((x, y) => sign(x.x - y.x) < 0)(intersections);
                
                // Disqualify multiple intersections
                foreach(j; 1..intersections.length)
                    if (intersections[j].x == intersections[j-1].x)
                        intersections[j].direction = intersections[j-1].direction = 0;
                
                // Inspect scanline and deduce orientations of intersected contours
                foreach(j; 0..intersections.length)
                    if (intersections[j].direction)
                        orientations[intersections[j].contourIndex] += 2*((j&1)^(intersections[j].direction > 0))-1;
                intersections.length = 0;
            }
        }
        // Reverse contours that have the opposite orientation
        foreach(i; 0..contours.length)
            if (orientations[i] < 0)
                contours[i].reverse();
    }
}