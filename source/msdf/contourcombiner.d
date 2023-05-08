module msdf.contourcombiner;

import std.algorithm : min, max;
import std.math : abs;

import inmath.linalg;
import msdf.contour;
import msdf.edgeselectors;
import msdf.shape;

private {
    void _msdfInitDistance(ref double distance) {
        distance = double.max;
    }

    void _msdfInitDistance(ref MultiDistance distance) {
        distance.r = double.max;
        distance.g = double.max;
        distance.b = double.max;
    }

    double _msdfResolveDistance(double distance) {
        return distance;
    }

    double _msdfResolveDistance( MultiDistance distance) {
        return max(min(distance.r, distance.g), min(max(distance.r, distance.g), distance.b));
    }
}

class SimpleContourCombiner(EdgeSelector) {
private:
    EdgeSelector shapeEdgeSelector;

public:
    alias EdgeSelectorType = EdgeSelector;
    alias DistanceType = EdgeSelector.DistanceType;

    this( Shape shape) { 
    }

    void reset( vec2d p) {
        shapeEdgeSelector.reset(p);
    }

    EdgeSelector edgeSelector(int i) {
        return shapeEdgeSelector;
    }

    DistanceType distance()  {
        return shapeEdgeSelector.distance();
    }
}

class OverlappingContourCombiner(EdgeSelector) {
private:
    vec2d p;
    int[] windings;
    EdgeSelector[] edgeSelectors;

public:
    alias EdgeSelectorType = EdgeSelector;
    alias DistanceType = EdgeSelector.DistanceType;

    this( Shape shape) {
        windings.reserve(shape.contours.length);
        foreach (contour; shape.contours)
            windings ~= contour.winding();
        edgeSelectors.length = shape.contours.length;
    }

    void reset( vec2d p) {
        this.p = p;

        foreach (EdgeSelector contourEdgeSelector; edgeSelectors)
            contourEdgeSelector.reset(p);
    }

    EdgeSelector edgeSelector(int i) {
        return edgeSelectors[i];
    }

    DistanceType distance()  {
        int contourCount = cast(int) edgeSelectors.length;

        EdgeSelector shapeEdgeSelector;
        EdgeSelector innerEdgeSelector;
        EdgeSelector outerEdgeSelector;
        shapeEdgeSelector.reset(p);
        innerEdgeSelector.reset(p);
        outerEdgeSelector.reset(p);
        
        for (int i = 0; i < contourCount; ++i) {
            DistanceType edgeDistance = edgeSelectors[i].distance();
            shapeEdgeSelector.merge(edgeSelectors[i]);
            if (windings[i] > 0 && _msdfResolveDistance(edgeDistance) >= 0)
                innerEdgeSelector.merge(edgeSelectors[i]);
            if (windings[i] < 0 && _msdfResolveDistance(edgeDistance) <= 0)
                outerEdgeSelector.merge(edgeSelectors[i]);
        }

        DistanceType shapeDistance = shapeEdgeSelector.distance();
        DistanceType innerDistance = innerEdgeSelector.distance();
        DistanceType outerDistance = outerEdgeSelector.distance();
        double innerScalarDistance = _msdfResolveDistance(innerDistance);
        double outerScalarDistance = _msdfResolveDistance(outerDistance);
        DistanceType distance;
        _msdfInitDistance(distance);

        int winding = 0;
        if (innerScalarDistance >= 0 && abs(innerScalarDistance) <= abs(outerScalarDistance)) {
            distance = innerDistance;
            winding = 1;
            for (int i = 0; i < contourCount; ++i)
                if (windings[i] > 0) {
                    DistanceType contourDistance = edgeSelectors[i].distance();
                    if (abs(_msdfResolveDistance(contourDistance)) < abs(outerScalarDistance) && _msdfResolveDistance(contourDistance) > _msdfResolveDistance(distance))
                        distance = contourDistance;
                }
        } else if (outerScalarDistance <= 0 && abs(outerScalarDistance) < abs(innerScalarDistance)) {
            distance = outerDistance;
            winding = -1;
            for (int i = 0; i < contourCount; ++i)
                if (windings[i] < 0) {
                    DistanceType contourDistance = edgeSelectors[i].distance();
                    if (abs(_msdfResolveDistance(contourDistance)) < abs(innerScalarDistance) && _msdfResolveDistance(contourDistance) < _msdfResolveDistance(distance))
                        distance = contourDistance;
                }
        } else
            return shapeDistance;

        for (int i = 0; i < contourCount; ++i)
            if (windings[i] != winding) {
                DistanceType contourDistance = edgeSelectors[i].distance();
                if (_msdfResolveDistance(contourDistance)*_msdfResolveDistance(distance) >= 0 && abs(_msdfResolveDistance(contourDistance)) < abs(_msdfResolveDistance(distance)))
                    distance = contourDistance;
            }
        if (_msdfResolveDistance(distance) == _msdfResolveDistance(shapeDistance))
            distance = shapeDistance;
        return distance;
    }
}