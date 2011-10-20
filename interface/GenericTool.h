/**
    @file GenericTool.h

    @brief Generic tools (C++ based) for analysis.

    @author Haryo Sumowidagdo <Suharyo.Sumowidagdo@cern.ch>

    @date Sat Aug 28 20:57:42 CEST 2010

    $Id: GenericTool.h,v 1.1 2011/01/20 15:26:15 haryo Exp $
*/


#ifndef wilken_GenericTool_h
#define wilken_GenericTool_h


namespace {

    /**
        A predicate which returns TRUE if the first object has
        greater \f$ p_{T} \f$ than the second object.

        @param T The object class/type.

        @param a The first object.

        @param b The second object.
     */
    template<class T>
    bool
    PtGreaterThan(const T& a, const T& b) {
        return a.pt() > b.pt() ;
    }

    /**
        A predicate which returns TRUE if the first object has
        less \f$ p_{T} \f$ than the second object.

        @param T The object class/type.

        @param a The first object.

        @param b The second object.
     */
    template<class T>
    bool
    PtLessThan(const T& a, const T& b) {
        return a.pt() < b.pt() ;
    }

    /**
        A predicate which returns TRUE if the first indexed object
        has an attribute which is greater than the second indexed object.

        @param T The attribute's class/type.

        @param a The first object.

        @param b The second object.
     */
    template<class T>
    bool
    IndexedQuantityGreaterThan(const std::pair<size_t,T> a, const std::pair<size_t,T> b) {
        return a.second > b.second ;
    }

    /**
        A predicate which returns TRUE if the first indexed object
        has an attribute which is less than the second indexed object.

        @param T The attribute's class/type.

        @param a The first object.

        @param b The second object.
    */
    template<class T>
    bool
    IndexedQuantityLessThan(const std::pair<size_t,T> a, const std::pair<size_t,T> b){
        return a.second < b.second ;
    }

    /**
        A predicate which returns TRUE if the first indexed object
        has an attribute which is greater in absolute value
        than the second indexed object.

        @param T The attribute's class/type.

        @param a The first object.

        @param b The second object.
     */
    template<class T>
    bool
    IndexedQuantityAbsGreaterThan(const std::pair<size_t,T> a, const std::pair<size_t,T> b) {
        return fabs(a.second) > fabs(b.second) ;
    }

    /**
        A predicate which returns TRUE if the first indexed object
        has an attribute which is less in absolute value
        than the second indexed object.

        @param T The attribute's class/type.

        @param a The first object.

        @param b The second object.
     */
    template<class T>
    bool
    IndexedQuantityAbsLessThan(const std::pair<size_t,T> a, const std::pair<size_t,T> b) {
        return fabs(a.second) < fabs(b.second) ;
    }

}

#endif
