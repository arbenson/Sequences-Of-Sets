const Set0 = NTuple{0,Int32}
const Set1 = NTuple{1,Int32}
const Set2 = NTuple{2,Int32}
const Set3 = NTuple{3,Int32}
const Set4 = NTuple{4,Int32}
const Set5 = NTuple{5,Int32}
const SmallFixedSizeSet = Union{Set0, Set1, Set2, Set3, Set4, Set5}

"""
SFSS_FromOrderedVec
----------

SFSSFromOrderedVec(v::Vector{Int32})

Create a small fixed-size set from a vector v

v must contain be of length at most 5
"""
function SFSSFromOrderedVec(v::Vector{Int32})
    if length(v) == 0; return Set0(); end
    if length(v) == 1; return Set1(v); end
    if length(v) == 2; return Set2(v); end
    if length(v) == 3; return Set3(v); end
    if length(v) == 4; return Set4(v); end
    if length(v) == 5; return Set5(v); end
    error("Set size is too large (maximum size is 5)")
end

SFSSFromOrderedVec(v::Vector{Int64}) =
    SFSSFromOrderedVec(convert(Vector{Int32}, v))
SFSSFromVec(v::Vector{Int64}) =
    SFSSFromOrderedVec(sort(v))
SFSSFromVec(v::Vector{Int32}) =
    SFSSFromOrderedVec(sort(v))
SFSSFromSet(S::Set{Int64}) =
    SFSSFromOrderedVec(sort(collect(S)))
SFSSFromSet(S::Set{Int32}) =
    SFSSFromOrderedVec(sort(collect(S)))

# Assumption: sets are sorted
function SFSS_issubset(A::SmallFixedSizeSet, B::SmallFixedSizeSet)
    nA, nB = length(A), length(B)
    Aind, Bind = 1, 1
    while Aind <= nA && Bind <= nB
        aval, bval = A[Aind], B[Bind]        
        if   aval < bval;  return false
        else Bind += 1;
        end
        Aind += (aval == bval)
    end
    return Aind > nA
end

"""
SFSS_intersect
-----------------

SFSS_intersect(A::SmallFixedSizeSet, B::SmallFixedSizeSet)

Returns the intersect of small fixed-size sets A and B as a small fixed-size set.
"""
function SFSS_intersect(A::SmallFixedSizeSet, B::SmallFixedSizeSet)
    intersect_vec = Int32[]
    Aind, Bind = 1, 1
    nA, nB = length(A), length(B)
    while Aind <= nA && Bind <= nB
        aval, bval = A[Aind], B[Bind]
        if aval == bval
            push!(intersect_vec, aval)
            Aind += 1
            Bind += 1
        elseif aval < bval; Aind += 1;
        else                Bind += 1;
        end
    end
    return SFSSFromOrderedVec(intersect_vec)
end

"""
SFSS_setdiff
------------

SFSS_setdiff(A::SmallFixedSizeSet, B::SmallFixedSizeSet)

Returns the elements in A that are not also in B as a small fixed-size set.
"""
function SFSS_setdiff(A::SmallFixedSizeSet, B::SmallFixedSizeSet)
    setdiff_vec = Int32[]
    Aind, Bind = 1, 1
    nA, nB = length(A), length(B)
    while Aind <= nA && Bind <= nB
        aval, bval = A[Aind], B[Bind]
        if aval == bval
            Aind += 1
            Bind += 1
        elseif aval < bval
            push!(setdiff_vec, aval)
            Aind += 1
        else
            Bind += 1
        end
    end
    while Aind <= nA
        push!(setdiff_vec, A[Aind])
        Aind += 1
    end
    return SFSSFromOrderedVec(setdiff_vec)
end

"""
SFSS_union
--------------

SFSS_union(A::SmallFixedSizeSet, B::SmallFixedSizeSet)

Returns the union of small fixed-size sets A and B.
"""
function SFSS_union(A::SmallFixedSizeSet, B::SmallFixedSizeSet)
    union_vec = Int32[]
    Aind, Bind = 1, 1
    nA, nB = length(A), length(B)
    while Aind <= nA && Bind <= nB
        aval, bval = A[Aind], B[Bind]
        if aval == bval
            push!(union_vec, aval)
            Aind += 1
            Bind += 1
        elseif aval < bval
            push!(union_vec, aval)
            Aind += 1
        else
            push!(union_vec, bval)
            Bind += 1
        end
    end
    while Aind <= nA
        push!(union_vec, A[Aind])
        Aind += 1
    end
    while Bind <= nB
        push!(union_vec, B[Bind])
        Bind += 1
    end
    return SFSSFromOrderedVec(union_vec)
end

