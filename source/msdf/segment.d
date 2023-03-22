module msdf.segment;
import inmath;
import inmath.interpolate;
import inmath.math;
import msdf.edgeholder;
import msdf.common;
import core.stdc.tgmath;

enum EdgeColor {
    BLACK,
    RED,
    GREEN,
    YELLOW,
    BLUE,
    MAGENTA,
    CYAN,
    WHITE
}

enum MSDFGEN_CUBIC_SEARCH_STARTS = 4.0;
enum MSDFGEN_CUBIC_SEARCH_STEPS = 4.0;

private {
    void _segPointBounds(vec2d p, ref double l, ref double b, ref double r, ref double t) {
        if (p.x < l) l = p.x;
        if (p.y < b) b = p.y;
        if (p.x > r) r = p.x;
        if (p.y > t) t = p.y;
    }
}

abstract
class EdgeSegment {
public:
    EdgeColor color;

    this(EdgeColor edgeColor = EdgeColor.WHITE) { color = edgeColor; }

    /// Creates a copy of the edge segment.
    abstract EdgeSegment clone() const; 
    /// Returns the point on the edge specified by the parameter (between 0 and 1).
    abstract vec2d point(double param) const;
    /// Returns the direction the edge has at the point specified by the parameter.
    abstract vec2d direction(double param) const;
    /// Returns the change of direction (second derivative) at the point specified by the parameter.
    abstract vec2d directionChange(double param) const;
    /// Returns the minimum signed distance between origin and the edge.
    abstract vec2d signedDistance(vec2d origin, ref double param) const;
    /// Outputs a list of (at most three) intersections (their X coordinates) with an infinite horizontal scanline at y and returns how many there are.
    abstract int scanlineIntersections(double[3] x, int[3] dy, double y) const;
    /// Adjusts the bounding box to fit the edge segment.
    abstract void bound(ref double l, ref double b, ref double r, ref double t) const;

    /// Reverses the edge (swaps its start point and end point).
    abstract void reverse();
    /// Moves the start point of the edge segment.
    abstract void moveStartPoint(vec2d to);
    /// Moves the end point of the edge segment.
    abstract void moveEndPoint(vec2d to);
    /// Splits the edge segments into thirds which together represent the original edge.
    abstract void splitInThirds(ref EdgeSegment part1, ref EdgeSegment part2, ref EdgeSegment part3) const;
}

class LinearSegment : EdgeSegment {
public:
    vec2d[2] p;

    this(vec2d p0, vec2d p1, EdgeColor edgeColor = EdgeColor.WHITE) {
        super(edgeColor);
        this.p[0] = p0;
        this.p[1] = p1;
        this.color = edgeColor;
    }

    override
    EdgeSegment clone() const {
        return new LinearSegment(p[0], p[1], color);
    }

    override
    vec2d point(double param) const {
        return lerp(p[0], p[1], param);
    }

    override
    vec2d direction(double param) const {
        return p[1]-p[0];
    }

    override
    vec2d directionChange(double param) const {
        return vec2d.zero;
    }

    override
    vec2d signedDistance(vec2d origin, ref double param) const {
        vec2d aq = origin-p[0];
        vec2d ab = p[1]-p[0];
        param = dot(aq, ab)/dot(ab, ab);
        vec2d eq = p[param > .5]-origin;
        double endpointDistance = eq.length();
        if (param > 0 && param < 1) {
            double orthoDistance = dot(ab.getOrthonormal(false), aq);
            if (abs(orthoDistance) < endpointDistance)
                return vec2d(orthoDistance, 0);
        }
        return vec2d(nzsign(cross(aq, ab))*endpointDistance, fabs(dot(ab.normalized(), eq.normalized())));
    }

    override
    int scanlineIntersections(double[3] x, int[3] dy, double y) const {
        if ((y >= p[0].y && y < p[1].y) || (y >= p[1].y && y < p[0].y)) {
            double param = (y-p[0].y)/(p[1].y-p[0].y);
            x[0] = lerp(p[0].x, p[1].x, param);
            dy[0] = sign(p[1].y-p[0].y);
            return 1;
        }
        return 0;
    }

    override
    void bound(ref double l, ref double b, ref double r, ref double t) const {
        _segPointBounds(p[0], l, b, r, t);
        _segPointBounds(p[1], l, b, r, t);
    }

    override
    void reverse() {
        vec2d tmp = p[0];
        p[0] = p[1];
        p[1] = tmp;
    }

    override
    void moveStartPoint(vec2d to) {
        p[0] = to;
    }

    override
    void moveEndPoint(vec2d to) {
        p[1] = to;
    }

    override 
    void splitInThirds(ref EdgeSegment part1, ref EdgeSegment part2, ref EdgeSegment part3) const {
        part = new LinearSegment(p[0], point(1/3.0), color);
        part = new LinearSegment(point(1/3.0), point(2/3.0), color);
        part = new LinearSegment(point(2/3.0), p[1], color);
    }

    vec2d length() const {
        return (p[0]-p[1]).length;
    }
}

class QuadraticSegment : EdgeSegment {
public:
    vec2d[3] p;

    this(vec2d p0, vec2d p1, vec2d p2, EdgeColor edgeColor = EdgeColor.WHITE) {
        super(edgeColor);

        if (p1 == p0 || p1 == p2)
            p1 = 0.5*(p0+p2);
        this.p[0] = p0;
        this.p[1] = p1;
        this.p[2] = p2;
        
    }

    override
    EdgeSegment clone() const {
        return new QuadraticSegment(p[0], p[1], p[2], color);
    }

    override
    vec2d point(double param) const {
        return lerp(lerp(p[0], p[1], param), lerp(p[1], p[2], param), param);
    }

    override
    vec2d direction(double param) const {
        vec2d tangent = lerp(p[1]-p[0], p[2]-p[1], param);
        if (!tangent)
            return p[2]-p[0];
        return tangent;
    }

    override
    vec2d directionChange(double param) const {
        return (p[2]-p[1])-(p[1]-p[0]);
    }

    override
    vec2d signedDistance(vec2d origin, ref double param) const {
        vec2d qa = p[0]-origin;
        vec2d ab = p[1]-p[0];
        vec2d br = p[2]-p[1]-ab;
        double a = dot(br, br);
        double b = 3*dot(ab, br);
        double c = 2*dot(ab, ab)+dot(qa, br);
        double d = dot(qa, ab);
        double[3] t;
        int solutions = solveCubic(t, a, b, c, d);

        vec2d epDir = direction(0);
        double minDistance = nzsign(cross(epDir, qa))*qa.length(); // distance from A
        param = -dot(qa, epDir)/dot(epDir, epDir);
        {
            epDir = direction(1);
            double distance = (p[2]-origin).length(); // distance from B
            if (distance < abs(minDistance)) {
                minDistance = nzsign(cross(epDir, p[2]-origin))*distance;
                param = dot(origin-p[1], epDir)/dot(epDir, epDir);
            }
        }
        for (int i = 0; i < solutions; ++i) {
            if (t[i] > 0 && t[i] < 1) {
                vec2d qe = qa+2*t[i]*ab+t[i]*t[i]*br;
                double distance = qe.length();
                if (distance <= abs(minDistance)) {
                    minDistance = nzsign(cross(ab+t[i]*br, qe))*distance;
                    param = t[i];
                }
            }
        }

        if (param >= 0 && param <= 1)
            return vec2d(minDistance, 0);
        if (param < .5)
            return vec2d(minDistance, abs(dot(direction(0).normalized(), qa.normalized())));
        else
            return vec2d(minDistance, abs(dot(direction(1).normalized(), (p[2]-origin).normalized())));
    }
    

    override
    int scanlineIntersections(double[3] x, int[3] dy, double y) const {
        int total = 0;
        int nextDY = y > p[0].y ? 1 : -1;
        x[total] = p[0].x;
        if (p[0].y == y) {
            if (p[0].y < p[1].y || (p[0].y == p[1].y && p[0].y < p[2].y))
                dy[total++] = 1;
            else
                nextDY = 1;
        }
        {
            Vector2 ab = p[1]-p[0];
            Vector2 br = p[2]-p[1]-ab;
            double t[2];
            int solutions = solveQuadratic(t, br.y, 2*ab.y, p[0].y-y);
            // Sort solutions
            double tmp;
            if (solutions >= 2 && t[0] > t[1])
                tmp = t[0], t[0] = t[1], t[1] = tmp;
            for (int i = 0; i < solutions && total < 2; ++i) {
                if (t[i] >= 0 && t[i] <= 1) {
                    x[total] = p[0].x+2*t[i]*ab.x+t[i]*t[i]*br.x;
                    if (nextDY*(ab.y+t[i]*br.y) >= 0) {
                        dy[total++] = nextDY;
                        nextDY = -nextDY;
                    }
                }
            }
        }
        if (p[2].y == y) {
            if (nextDY > 0 && total > 0) {
                --total;
                nextDY = -1;
            }
            if ((p[2].y < p[1].y || (p[2].y == p[1].y && p[2].y < p[0].y)) && total < 2) {
                x[total] = p[2].x;
                if (nextDY < 0) {
                    dy[total++] = -1;
                    nextDY = 1;
                }
            }
        }
        if (nextDY != (y >= p[2].y ? 1 : -1)) {
            if (total > 0)
                --total;
            else {
                if (abs(p[2].y-y) < abs(p[0].y-y))
                    x[total] = p[2].x;
                dy[total++] = nextDY;
            }
        }
        return total;
    }

    override
    void bound(ref double l, ref double b, ref double r, ref double t) const {
        _segPointBounds(p[0], l, b, r, t);
        _segPointBounds(p[2], l, b, r, t);
        vec2d bot = (p[1]-p[0])-(p[2]-p[1]);
        if (bot.x) {
            double param = (p[1].x-p[0].x)/bot.x;
            if (param > 0 && param < 1)
                _segPointBounds(point(param), l, b, r, t);
        }
        if (bot.y) {
            double param = (p[1].y-p[0].y)/bot.y;
            if (param > 0 && param < 1)
                _segPointBounds(point(param), l, b, r, t);
        }
    }
    
    override
    void reverse() {
        vec2d tmp = p[0];
        p[0] = p[2];
        p[2] = tmp;
    }

    override
    void moveStartPoint(vec2d to) {
        p[0] = to;
    }

    override
    void moveEndPoint(vec2d to) {
        vec2d origEDir = p[2]-p[1];
        vec2d origP1 = p[1];
        p[1] += cross(p[2]-p[1], to-p[2])/cross(p[2]-p[1], p[0]-p[1])*(p[0]-p[1]);
        p[2] = to;
        if (dot(origEDir, p[2]-p[1]) < 0)
            p[1] = origP1;
    }

    double length() const {
        vec2d ab = p[1]-p[0];
        vec2d br = p[2]-p[1]-ab;
        double abab = dot(ab, ab);
        double abbr = dot(ab, br);
        double brbr = dot(br, br);
        double abLen = sqrt(abab);
        double brLen = sqrt(brbr);
        double crs = cross(ab, br);
        double h = sqrt(abab+abbr+abbr+brbr);
        return (
            brLen*((abbr+brbr)*h-abbr*abLen)+
            crs*crs*log((brLen*h+abbr+brbr)/(brLen*abLen+abbr))
        )/(brbr*brLen);
    }

    EdgeSegment convertToCubic() const {
        return new CubicSegment(p[0], lerp(p[0], p[1], 2/3.), lerp(p[1], p[2], 1/3.), p[2], color);
    }
    
}

class CubicSegment : EdgeSegment {
public:
    vec2d[4] p;

    this(vec2d p0, vec2d p1, vec2d p2, vec2d p3, EdgeColor edgeColor = EdgeColor.WHITE) {
        super(edgeColor);

        if ((p1 == p0 || p1 == p3) && (p2 == p0 || p2 == p3)) {
            p1 = lerp(p0, p3, 1/3.);
            p2 = lerp(p0, p3, 2/3.);
        }
        this.p[0] = p0;
        this.p[1] = p1;
        this.p[2] = p2;
        this.p[3] = p3;
    }

    override
    EdgeSegment clone() const {
        return new CubicSegment(p[0], p[1], p[2], p[3], color);
    }

    override
    vec2d point(double param) const {
        vec2d p12 = lerp(p[1], p[2], param);
        return lerp(lerp(lerp(p[0], p[1], param), p12, param), lerp(p12, lerp(p[2], p[3], param), param), param);
    }

    override
    vec2d direction(double param) const {
        vec2d tangent = lerp(lerp(p[1]-p[0], p[2]-p[1], param), lerp(p[2]-p[1], p[3]-p[2], param), param);
        if (!tangent) {
            if (param == 0) return p[2]-p[0];
            if (param == 1) return p[3]-p[1];
        }
        return tangent;
    }

    override
    vec2d directionChange(double param) const {
        return lerp((p[2]-p[1])-(p[1]-p[0]), (p[3]-p[2])-(p[2]-p[1]), param);
    }

    override
    vec2d signedDistance(vec2d origin, ref double param) const {
        vec2d qa = p[0]-origin;
        vec2d ab = p[1]-p[0];
        vec2d br = p[2]-p[1]-ab;
        vec2d as = (p[3]-p[2])-(p[2]-p[1])-br;

        vec2d epDir = direction(0);
        double minDistance = nzsign(cross(epDir, qa))*qa.length(); // distance from A
        param = -dot(qa, epDir)/dot(epDir, epDir);
        {
            epDir = direction(1);
            double distance = (p[3]-origin).length(); // distance from B
            if (distance < abs(minDistance)) {
                minDistance = nzsign(cross(epDir, p[3]-origin))*distance;
                param = dot(epDir-(p[3]-origin), epDir)/dot(epDir, epDir);
            }
        }
        // Iterative minimum distance search
        for (int i = 0; i <= MSDFGEN_CUBIC_SEARCH_STARTS; ++i) {
            double t = (double) i/MSDFGEN_CUBIC_SEARCH_STARTS;
            vec2d qe = qa+3*t*ab+3*t*t*br+t*t*t*as;
            for (int step = 0; step < MSDFGEN_CUBIC_SEARCH_STEPS; ++step) {
                // Improve t
                vec2d d1 = 3*ab+6*t*br+3*t*t*as;
                vec2d d2 = 6*br+6*t*as;
                t -= dot(qe, d1)/(dot(d1, d1)+dot(qe, d2));
                if (t <= 0 || t >= 1)
                    break;
                qe = qa+3*t*ab+3*t*t*br+t*t*t*as;
                double distance = qe.length();
                if (distance < abs(minDistance)) {
                    minDistance = nzsign(cross(d1, qe))*distance;
                    param = t;
                }
            }
        }

        if (param >= 0 && param <= 1)
            return vec2d(minDistance, 0);
        if (param < .5)
            return vec2d(minDistance, abs(dot(direction(0).normalized(), qa.normalized())));
        else
            return vec2d(minDistance, abs(dot(direction(1).normalized(), (p[3]-origin).normalized())));
    }

    override
    int scanlineIntersections(double[3] x, int[3] dy, double y) const {
        int total = 0;
        int nextDY = y > p[0].y ? 1 : -1;
        x[total] = p[0].x;
        if (p[0].y == y) {
            if (p[0].y < p[1].y || (p[0].y == p[1].y && (p[0].y < p[2].y || (p[0].y == p[2].y && p[0].y < p[3].y))))
                dy[total++] = 1;
            else
                nextDY = 1;
        }
        {
            vec2d ab = p[1]-p[0];
            vec2d br = p[2]-p[1]-ab;
            vec2d as = (p[3]-p[2])-(p[2]-p[1])-br;
            double t[3];
            int solutions = solveCubic(t, as.y, 3*br.y, 3*ab.y, p[0].y-y);
            // Sort solutions
            double tmp;
            if (solutions >= 2) {
                if (t[0] > t[1])
                    tmp = t[0], t[0] = t[1], t[1] = tmp;
                if (solutions >= 3 && t[1] > t[2]) {
                    tmp = t[1], t[1] = t[2], t[2] = tmp;
                    if (t[0] > t[1])
                        tmp = t[0], t[0] = t[1], t[1] = tmp;
                }
            }
            for (int i = 0; i < solutions && total < 3; ++i) {
                if (t[i] >= 0 && t[i] <= 1) {
                    x[total] = p[0].x+3*t[i]*ab.x+3*t[i]*t[i]*br.x+t[i]*t[i]*t[i]*as.x;
                    if (nextDY*(ab.y+2*t[i]*br.y+t[i]*t[i]*as.y) >= 0) {
                        dy[total++] = nextDY;
                        nextDY = -nextDY;
                    }
                }
            }
        }
        if (p[3].y == y) {
            if (nextDY > 0 && total > 0) {
                --total;
                nextDY = -1;
            }
            if ((p[3].y < p[2].y || (p[3].y == p[2].y && (p[3].y < p[1].y || (p[3].y == p[1].y && p[3].y < p[0].y)))) && total < 3) {
                x[total] = p[3].x;
                if (nextDY < 0) {
                    dy[total++] = -1;
                    nextDY = 1;
                }
            }
        }
        if (nextDY != (y >= p[3].y ? 1 : -1)) {
            if (total > 0)
                --total;
            else {
                if (abs(p[3].y-y) < abs(p[0].y-y))
                    x[total] = p[3].x;
                dy[total++] = nextDY;
            }
        }
        return total;
    }

    override
    void bound(ref double l, ref double b, ref double r, ref double t) const {
        _segPointBounds(p[0], l, b, r, t);
        _segPointBounds(p[3], l, b, r, t);
        vec2d a0 = p[1]-p[0];
        vec2d a1 = 2*(p[2]-p[1]-a0);
        vec2d a2 = p[3]-3*p[2]+3*p[1]-p[0];
        double[2] params;
        int solutions;
        solutions = solveQuadratic(params, a2.x, a1.x, a0.x);
        for (int i = 0; i < solutions; ++i)
            if (params[i] > 0 && params[i] < 1)
                _segPointBounds(point(params[i]), l, b, r, t);
        solutions = solveQuadratic(params, a2.y, a1.y, a0.y);
        for (int i = 0; i < solutions; ++i)
            if (params[i] > 0 && params[i] < 1)
                _segPointBounds(point(params[i]), l, b, r, t);
    }

    override
    void reverse() {
        vec2d tmp = p[0];
        p[0] = p[3];
        p[3] = tmp;
        tmp = p[1];
        p[1] = p[2];
        p[2] = tmp;
    }

    override
    void moveStartPoint(vec2d to) {
        p[1] += to-p[0];
        p[0] = to;
    }

    override
    void moveEndPoint(vec2d to) {
        p[2] += to-p[3];
        p[3] = to;
    }

    override
    void splitInThirds(ref EdgeSegment part1, ref EdgeSegment part2, ref EdgeSegment part3) const {
        part1 = new CubicSegment(p[0], p[0] == p[1] ? p[0] : lerp(p[0], p[1], 1/3.), lerp(lerp(p[0], p[1], 1/3.), lerp(p[1], p[2], 1/3.), 1/3.), point(1/3.), color);
        part2 = new CubicSegment(point(1/3.),
            lerp(lerp(lerp(p[0], p[1], 1/3.), lerp(p[1], p[2], 1/3.), 1/3.), lerp(lerp(p[1], p[2], 1/3.), lerp(p[2], p[3], 1/3.), 1/3.), 2/3.),
            lerp(lerp(lerp(p[0], p[1], 2/3.), lerp(p[1], p[2], 2/3.), 2/3.), lerp(lerp(p[1], p[2], 2/3.), lerp(p[2], p[3], 2/3.), 2/3.), 1/3.),
            point(2/3.), color);
        part3 = new CubicSegment(point(2/3.), lerp(lerp(p[1], p[2], 2/3.), lerp(p[2], p[3], 2/3.), 2/3.), p[2] == p[3] ? p[3] : lerp(p[2], p[3], 2/3.), p[3], color);
    }

    void deconverge(int param, double amount) {
        vec2d dir = direction(param);
        vec2d normal = dir.getOrthonormal();
        double h = dot(directionChange(param)-dir, normal);
        switch (param) {
            case 0:
                p[1] += amount*(dir+sign(h)*sqrt(abs(h))*normal);
                break;
            case 1:
                p[2] -= amount*(dir-sign(h)*sqrt(abs(h))*normal);
                break;
        }
    }
}