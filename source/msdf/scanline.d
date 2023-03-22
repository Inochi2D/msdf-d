module msdf.scanline;
import inmath;
import std.algorithm.sorting;

/// Fill rule dictates how intersection total is interpreted during rasterization.
enum FillRule {
    FILL_NONZERO,
    FILL_ODD, // even-odd
    FILL_POSITIVE,
    FILL_NEGATIVE
}

private {
    double _msdfOverlap(Scanline a, Scanline b, double xFrom, double xTo, FillRule fillRule) {
        double total = 0;
        bool aInside = false, bInside = false;
        int ai = 0, bi = 0;
        double ax = a.intersections.length != 0? a.intersections[ai].x : xTo;
        double bx = b.intersections.length != 0 ? b.intersections[bi].x : xTo;
        while (ax < xFrom || bx < xFrom) {
            double xNext = min(ax, bx);
            if (ax == xNext && ai < cast(int)a.intersections.length) {
                aInside = _msdfInterpretFillRule(a.intersections[ai].direction, fillRule);
                ax = ++ai < cast(int)a.intersections.length ? a.intersections[ai].x : xTo;
            }
            if (bx == xNext && bi < cast(int)b.intersections.length) {
                bInside = _msdfInterpretFillRule(b.intersections[bi].direction, fillRule);
                bx = ++bi < cast(int)b.intersections.length ? b.intersections[bi].x : xTo;
            }
        }
        double x = xFrom;
        while (ax < xTo || bx < xTo) {
            double xNext = min(ax, bx);
            if (aInside == bInside)
                total += xNext-x;
            if (ax == xNext && ai < cast(int)a.intersections.length) {
                aInside = _msdfInterpretFillRule(a.intersections[ai].direction, fillRule);
                ax = ++ai < cast(int)a.intersections.length ? a.intersections[ai].x : xTo;
            }
            if (bx == xNext && bi < cast(int)b.intersections.length) {
                bInside = _msdfInterpretFillRule(b.intersections[bi].direction, fillRule);
                bx = ++bi < cast(int)b.intersections.length ? b.intersections[bi].x : xTo;
            }
            x = xNext;
        }
        if (aInside == bInside)
            total += xTo-x;
        return total;
    }
    
    bool _msdfInterpretFillRule(int intersections, FillRule fillRule) {
        switch (fillRule) {
            case FillRule.FILL_NONZERO:
                return intersections != 0;
            case FillRule.FILL_ODD:
                return intersections&1;
            case FillRule.FILL_POSITIVE:
                return intersections > 0;
            case FillRule.FILL_NEGATIVE:
                return intersections < 0;
            default: return false;
        }
    }
}

/// Represents a horizontal scanline intersecting a shape.
class Scanline {
private:
    Intersection[] intersections;
    int lastIndex;

    void preprocess() {
        lastIndex = 0;
        if (intersections.length > 0) {
            sort!((x, y) => sign(x.x - y.x) < 0)(intersections);
            int totalDirection = 0;
            foreach(ref intersection; intersections) {
                totalDirection += intersection.direction;
                intersection.direction = totalDirection;
            }
        }
    }

    int moveTo(double x) {
        if (intersections.length == 0) return -1;

        int index = lastIndex;
        if (x < intersections[index].x) {
            do {
                if (index == 0) {
                    lastIndex = 0;
                    return -1;
                }
                --index;
            } while (x < intersections[index].x);
        } else {
            while (index < cast(int)intersections.length-1 && x >= intersections[index+1].x)
                ++index;
        }
        lastIndex = index;
        return index;
    }

public:

    /// An intersection with the scanline.
    struct Intersection {

        /// X coordinate.
        double x;

        /// Normalized Y direction of the oriented edge at the point of intersection.
        int direction;
    }

    this() {
        lastIndex = 0;
    }

    /// Populates the intersection list.
    void setIntersections(Intersection[] intersections) {
        this.intersections = intersections;
        this.preprocess();
    }

    /// Returns the number of intersections left of x.
    int countIntersection(double x) {
        return moveTo(x)+1;
    }

    /// Returns the total sign of intersections left of x.
    int sumIntersections(double x) {
        int index = moveTo(x);
        if (index >= 0)
            return intersections[index].direction;
        return 0;
    }

    /// Decides whether the scanline is filled at x based on fill rule.
    bool filled(double x, FillRule fillRule) {
        return _msdfInterpretFillRule(sumIntersections(x), fillRule);
    }
}