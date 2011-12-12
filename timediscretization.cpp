#include "timediscretization.hpp"


TimeDiscretization::TimeDiscretization(World& _W, Potential& _Pot, ObserverXYZ &_O) : W(_W), Pot(_Pot), O(_O)
{
    // TimeDiscretization initializer list initializes directly in the order of their declaration in .hpp
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
