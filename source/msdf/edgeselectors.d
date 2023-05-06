module msdf.edgeselectors;

import std.typecons : Rebindable;
import std.math : abs;
import inmath.linalg;
import msdf.signeddistance;
import msdf.segment;

enum DISTANCE_DELTA_FACTOR = 1.001;

struct MultiDistance {
    double r, g, b;
}

class PseudoDistanceSelectorBase {
private:
    SignedDistance minTrueDistance;
    double minNegativePseudoDistance;
    double minPositivePseudoDistance;
    Rebindable!(const EdgeSegment) nearEdge;
    double nearEdgeParam;

public:
    struct EdgeCache {
        vec2d point;
        double absDistance = 0;
        double aDomainDistance = 0, bDomainDistance = 0;
        double aPseudoDistance = 0, bPseudoDistance = 0;
    }

    static bool getPseudoDistance(ref double distance, in vec2d ep, in vec2d edgeDir) {
        double ts = ep.dot(edgeDir);

        if (ts > 0) {
            double pseudoDistance = ep.cross(edgeDir);
            if (abs(pseudoDistance) < abs(distance)) {
                distance = pseudoDistance;
                return true;
            }
        }

        return false;
    }

    this() {
        minNegativePseudoDistance = -abs(minTrueDistance.distance);
        minPositivePseudoDistance = abs(minTrueDistance.distance);
        nearEdge = null;
        nearEdgeParam = 0;
    }

    void reset(double delta) {
        minTrueDistance += (minTrueDistance.distance < 0 ? -1 : 1) * delta;
        minNegativePseudoDistance = -abs(minTrueDistance.distance);
        minPositivePseudoDistance = abs(minTrueDistance.distance);
        nearEdge = null;
        nearEdgeParam = 0;
    }

    bool isEdgeRelevant(in EdgeCache cache, const(EdgeSegment) edge, in vec2d p) const {
        double delta = DISTANCE_DELTA_FACTOR * (p-cache.point).length();

        return (
            cache.absDistance-delta <= abs(minTrueDistance.distance) ||
            abs(cache.aDomainDistance) < delta ||
            abs(cache.bDomainDistance) < delta ||
            (cache.aDomainDistance > 0 && (cache.aPseudoDistance < 0 ?
                cache.aPseudoDistance+delta >= minNegativePseudoDistance :
                cache.aPseudoDistance-delta <= minPositivePseudoDistance
            )) ||
            (cache.bDomainDistance > 0 && (cache.bPseudoDistance < 0 ?
                cache.bPseudoDistance+delta >= minNegativePseudoDistance :
                cache.bPseudoDistance-delta <= minPositivePseudoDistance
            ))
        );
    }

    void addEdgeTrueDistance(const(EdgeSegment) edge, in SignedDistance distance, double param) {
        if (distance < minTrueDistance) {
            minTrueDistance = distance;
            nearEdge = edge;
            nearEdgeParam = param;
        }
    }

    void addEdgePseudoDistance(double distance) {
        if (distance <= 0 && distance > minNegativePseudoDistance)
            minNegativePseudoDistance = distance;
        if (distance >= 0 && distance < minPositivePseudoDistance)
            minPositivePseudoDistance = distance;
    }

    void merge(in PseudoDistanceSelectorBase other) {
        if (other.minTrueDistance < minTrueDistance) {
            minTrueDistance = other.minTrueDistance;
            nearEdge = other.nearEdge;
            nearEdgeParam = other.nearEdgeParam;
        }
        if (other.minNegativePseudoDistance > minNegativePseudoDistance)
            minNegativePseudoDistance = other.minNegativePseudoDistance;
        if (other.minPositivePseudoDistance < minPositivePseudoDistance)
            minPositivePseudoDistance = other.minPositivePseudoDistance;
    }

    double computeDistance(in vec2d p) const {
        double minDistance = minTrueDistance.distance < 0 ? minNegativePseudoDistance : minPositivePseudoDistance;
        if (nearEdge) {
            SignedDistance distance = minTrueDistance;
            nearEdge.distanceToPseudoDistance(distance, p, nearEdgeParam);
            if (abs(distance.distance) < abs(minDistance))
                minDistance = distance.distance;
        }
        return minDistance;
    }

    SignedDistance trueDistance() const {
        return minTrueDistance;
    }
}

class PseudoDistanceSelector : PseudoDistanceSelectorBase { 
private:
    vec2d p;

public:
    alias DistanceType = double;

    void reset(in vec2d p) {
        double delta = DISTANCE_DELTA_FACTOR*(p-this.p).length();
        super.reset(delta);
        this.p = p;
    }

    void addEdge(ref EdgeCache cache, const(EdgeSegment) prevEdge, const(EdgeSegment) edge, const(EdgeSegment) nextEdge) {
        if (isEdgeRelevant(cache, edge, p)) {
            double param;
            SignedDistance distance = edge.signedDistance(p, param);
            addEdgeTrueDistance(edge, distance, param);
            cache.point = p;
            cache.absDistance = abs(distance.distance);

            vec2d ap = p-edge.point(0);
            vec2d bp = p-edge.point(1);
            vec2d aDir = edge.direction(0).normalize();
            vec2d bDir = edge.direction(1).normalize();
            vec2d prevDir = prevEdge.direction(1).normalize();
            vec2d nextDir = nextEdge.direction(0).normalize();
            double add = ap.dot((prevDir+aDir).normalize());
            double bdd = -bp.dot((bDir+nextDir).normalize());
            if (add > 0) {
                double pd = distance.distance;
                if (getPseudoDistance(pd, ap, -aDir))
                    addEdgePseudoDistance(pd = -pd);
                cache.aPseudoDistance = pd;
            }
            if (bdd > 0) {
                double pd = distance.distance;
                if (getPseudoDistance(pd, bp, bDir))
                    addEdgePseudoDistance(pd);
                cache.bPseudoDistance = pd;
            }
            cache.aDomainDistance = add;
            cache.bDomainDistance = bdd;
        }
    }

    DistanceType distance() const {
        return computeDistance(p);
    }
}

class MultiDistanceSelector {
private:
    vec2d p;
    PseudoDistanceSelectorBase r, g, b;

public:
    alias DistanceType = MultiDistance;
    alias EdgeCache = PseudoDistanceSelectorBase.EdgeCache;

    void reset(in vec2d p) {
        double delta = DISTANCE_DELTA_FACTOR*(p-this.p).length();
        r.reset(delta);
        g.reset(delta);
        b.reset(delta);
        this.p = p;
    }

    void addEdge(ref EdgeCache cache, const(EdgeSegment) prevEdge, const(EdgeSegment) edge, const(EdgeSegment) nextEdge) {
        if (
            (edge.color&EdgeColor.RED && r.isEdgeRelevant(cache, edge, p)) ||
            (edge.color&EdgeColor.GREEN && g.isEdgeRelevant(cache, edge, p)) ||
            (edge.color&EdgeColor.BLUE && b.isEdgeRelevant(cache, edge, p))
        ) {
            double param;
            SignedDistance distance = edge.signedDistance(p, param);
            if (edge.color&EdgeColor.RED)
                r.addEdgeTrueDistance(edge, distance, param);
            if (edge.color&EdgeColor.GREEN)
                g.addEdgeTrueDistance(edge, distance, param);
            if (edge.color&EdgeColor.BLUE)
                b.addEdgeTrueDistance(edge, distance, param);
            cache.point = p;
            cache.absDistance = fabs(distance.distance);

            vec2d ap = p-edge.point(0);
            vec2d bp = p-edge.point(1);
            vec2d aDir = edge.direction(0).normalize();
            vec2d bDir = edge.direction(1).normalize();
            vec2d prevDir = prevEdge.direction(1).normalize();
            vec2d nextDir = nextEdge.direction(0).normalize();
            double add = ap.dot((prevDir+aDir).normalize());
            double bdd = -bp.dot((bDir+nextDir).normalize());
            if (add > 0) {
                double pd = distance.distance;
                if (PseudoDistanceSelectorBase.getPseudoDistance(pd, ap, -aDir)) {
                    pd = -pd;
                    if (edge.color&EdgeColor.RED)
                        r.addEdgePseudoDistance(pd);
                    if (edge.color&EdgeColor.GREEN)
                        g.addEdgePseudoDistance(pd);
                    if (edge.color&EdgeColor.BLUE)
                        b.addEdgePseudoDistance(pd);
                }
                cache.aPseudoDistance = pd;
            }
            if (bdd > 0) {
                double pd = distance.distance;
                if (PseudoDistanceSelectorBase.getPseudoDistance(pd, bp, bDir)) {
                    if (edge.color&EdgeColor.RED)
                        r.addEdgePseudoDistance(pd);
                    if (edge.color&EdgeColor.GREEN)
                        g.addEdgePseudoDistance(pd);
                    if (edge.color&EdgeColor.BLUE)
                        b.addEdgePseudoDistance(pd);
                }
                cache.bPseudoDistance = pd;
            }
            cache.aDomainDistance = add;
            cache.bDomainDistance = bdd;
        }
    }

    void merge(in MultiDistanceSelector other) {
        r.merge(other.r);
        g.merge(other.g);
        b.merge(other.b);
    }

    DistanceType distance() const {
        MultiDistance multiDistance;
        multiDistance.r = r.computeDistance(p);
        multiDistance.g = g.computeDistance(p);
        multiDistance.b = b.computeDistance(p);
        return multiDistance;
    }

    SignedDistance trueDistance() const {
        SignedDistance distance = r.trueDistance();
        if (g.trueDistance() < distance)
            distance = g.trueDistance();
        if (b.trueDistance() < distance)
            distance = b.trueDistance();
        return distance;
    }
}