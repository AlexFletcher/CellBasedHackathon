/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "LongRangeSignallingExamplePde.hpp"
#include "Exception.hpp"

template<unsigned DIM>
LongRangeSignallingExamplePde<DIM>::LongRangeSignallingExamplePde(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
                                                            double duDtCoefficient,
                                                            double diffusionCoefficient,
                                                            double sourceCoefficient,
                                                            double sinkCoefficient)
    : mrCellPopulation(rCellPopulation),
      mDuDtCoefficient(duDtCoefficient),
      mDiffusionCoefficient(diffusionCoefficient),
      mSourceCoefficient(sourceCoefficient),
      mSinkCoefficient(sinkCoefficient)
{
}

template<unsigned DIM>
const AbstractCellPopulation<DIM,DIM>& LongRangeSignallingExamplePde<DIM>::rGetCellPopulation() const
{
    return mrCellPopulation;
}

template<unsigned DIM>
double LongRangeSignallingExamplePde<DIM>::ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& )
{
    return mDuDtCoefficient;
}

template<unsigned DIM>
double LongRangeSignallingExamplePde<DIM>::ComputeSourceTerm(const ChastePoint<DIM>& rX, double u, Element<DIM,DIM>* pElement)
{
    NEVER_REACHED;
    return 0.0;
}

template<unsigned DIM>
double LongRangeSignallingExamplePde<DIM>::ComputeSourceTermAtNode(const Node<DIM>& rNode, double u)
{
    double source_term = 0.0;

    // here do: if mrCellPopulation.GetCellCorrespondingToNode(rNode).HasLabel do sink and source
    // else: do sink only
    // Todo: define which property distinguishes our red and green cells (I would guess we would use
    // cell labels.
    std::set<unsigned> neighbour_containing_elements =
            GetNode(neighbour_location_index)->rGetContainingElementIndices();

    unsigned elem_index = *(element_indices.begin());
    CellPtr p_cell = this->GetCellUsingLocationIndex(elem_index);
    if (cell_iter->HasCellProperty<CellLabel>())
    {
        source_term = mSourceCoefficient - mSinkCoefficient*u;
    }
    else
    {
        source_term = -mSinkCoefficient*u;
    }

    return source_term;
}

template<unsigned DIM>
c_matrix<double,DIM,DIM> LongRangeSignallingExamplePde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    return mDiffusionCoefficient*identity_matrix<double>(DIM);
}

// Explicit instantiation
template class LongRangeSignallingExamplePde<1>;
template class LongRangeSignallingExamplePde<2>;
template class LongRangeSignallingExamplePde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LongRangeSignallingExamplePde)
