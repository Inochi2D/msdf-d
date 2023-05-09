module msdf.contour;
import msdf.segment;
import msdf.common;
import inmath;

private {
    void _msdfBoundPoint(ref double l, ref double b, ref double r, ref double t, vec2d p) {
        if (p.x < l) l = p.x;
        if (p.y < b) b = p.y;
        if (p.x > r) r = p.x;
        if (p.y > t) t = p.y;
    }

    double _msdfShoelace(vec2d a, vec2d b) {
        return (b.x-a.x)*(a.y+b.y);
    }
}

struct Contour {
public:
    EdgeSegment[] edges;

    void addEdge(EdgeSegment segment) nothrow {
        this.edges ~= segment;
    }

    void bound(ref double l, ref double b, ref double r, ref double t)  {
        foreach(edge; edges) {
            edge.bound(l, b, r, t);
        }
    }

    void boundMiters(ref double l, ref double b, ref double r, ref double t, double border, double miterLimit, int polarity)  {
        if (edges.length == 0) return;

        vec2d prevDir = edges[$-1].direction(1).normalized;
        foreach(edge; edges) {
            vec2d dir = -(edge.direction(0).normalized());
            if (polarity*cross(prevDir, dir) >= 0) {
                double miterLength = miterLimit;
                double q = .5*(1-dot(prevDir, dir));
                if (q > 0)
                    miterLength = min(1/sqrt(q), miterLimit);
                vec2d miter = edge.point(0)+border*miterLength*(prevDir+dir).normalized;
                _msdfBoundPoint(l, b, r, t, miter);
            }
            prevDir = edge.direction(1).normalized;
        }
    }

    int winding()  {
        if (edges.length == 0)
            return 0;
        double total = 0;
        if (edges.length == 1) {
            vec2d a = edges[0].point(0), b = edges[0].point(1/3.), c = edges[0].point(2/3.);
            total += _msdfShoelace(a, b);
            total += _msdfShoelace(b, c);
            total += _msdfShoelace(c, a);
        } else if (edges.length == 2) {
            vec2d a = edges[0].point(0), b = edges[0].point(.5), c = edges[1].point(0), d = edges[1].point(.5);
            total += _msdfShoelace(a, b);
            total += _msdfShoelace(b, c);
            total += _msdfShoelace(c, d);
            total += _msdfShoelace(d, a);
        } else {
            vec2d prev = edges[$-1].point(0);
            foreach (edge; edges) {
                vec2d cur = edges[$-1].point(0);
                total += _msdfShoelace(prev, cur);
                prev = cur;
            }
        }
        return cast(int)sign(total);
    }

    void reverse() {
        EdgeSegment tmp;
        for(int i = cast(int)edges.length/2; i > 0; --i) {
            tmp = edges[i-1];
            edges[i-1] = edges[edges.length-i];
            edges[edges.length-i] = tmp;
        }

        foreach(ref edge; edges) {
            edge.reverse();
        }
    }
}