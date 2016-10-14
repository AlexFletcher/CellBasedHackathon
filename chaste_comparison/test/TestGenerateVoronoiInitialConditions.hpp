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

#ifndef TESTGENERATEVORONOIINITIALCONDITIONS_HPP_
#define TESTGENERATEVORONOIINITIALCONDITIONS_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath>

#include "VoronoiVertexMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"
#include "VertexMeshWriter.hpp"
#include "VertexMeshReader.hpp"
#include "PottsMesh.hpp"
#include "PottsMeshWriter.hpp"
#include "PottsMeshGenerator.hpp"
#include "Debug.hpp"

#include "FakePetscSetup.hpp"


void GeneratePottsFromVertexMesh(PottsMesh<2>* pPottsMesh, VertexMesh<2,2>*  pVertexMesh)
{
   
    PRINT_2_VARIABLES(pPottsMesh->GetNumNodes(),pPottsMesh->GetNumElements());
 
    unsigned element_counter = 0;
    for (VertexMesh<2,2>::VertexElementIterator iter = pVertexMesh->GetElementIteratorBegin();
         iter != pVertexMesh->GetElementIteratorEnd();
         ++iter)
    {
        unsigned element_index = iter->GetIndex();

        std::vector<Node<2>*> nodes_in_element;
        
        for (PottsMesh<2>::NodeIterator node_iter = pPottsMesh->GetNodeIteratorBegin();
             node_iter != pPottsMesh->GetNodeIteratorEnd();
             ++node_iter)
        {
            unsigned node_index = node_iter->GetIndex();

            // you need to make the ElementIncludesPoint method public for this to work.
            if (pVertexMesh->ElementIncludesPoint(node_iter->rGetLocation(),element_index))
            {
                nodes_in_element.push_back(pPottsMesh->GetNode(node_index));
            }           
        }

        assert(nodes_in_element.size()>0);

        pPottsMesh->AddElement(new PottsElement<2>(element_counter, nodes_in_element));

        element_counter++;
    }

    PRINT_2_VARIABLES(pPottsMesh->GetNumNodes(),pPottsMesh->GetNumElements());
}


bool CustomComparisonForPair(std::pair<unsigned, double> pair_a, std::pair<unsigned, double> pair_b)
{
    return (pair_a.second > pair_b.second);
}


class TestGenerateVoronoiInitialConditions : public CxxTest::TestSuite
{
public:

    void TestGenerateSquareVoronoi() throw(Exception)
    {
#if BOOST_VERSION >= 105200

        // Generate a mesh that is 10 cells wide, 10 high, with 3 Lloyd's relaxation steps and target average element area 0.5
        VoronoiVertexMeshGenerator generator(10, 10, 10, 0.5);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMeshAfterReMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 100u);

        // Check average cell area is correct
        double average_area = 0.0;
        for (unsigned elem_idx = 0 ; elem_idx < p_mesh->GetNumElements() ; elem_idx++)
        {
            average_area += p_mesh->GetVolumeOfElement(elem_idx);
        }
        average_area /= double(p_mesh->GetNumElements());

        TS_ASSERT_DELTA(average_area, 0.5, 1e-6);

        // Create a vertex mesh writer
        VertexMeshWriter<2,2> vertex_mesh_writer("TestGenerateSquareVoronoi", "voronoi_square_100");
        vertex_mesh_writer.WriteFilesUsingMesh(*p_mesh);
        vertex_mesh_writer.WriteVtkUsingMesh(*p_mesh,"hi");

        c_vector<double, 2> upper_corner = p_mesh->CalculateBoundingBox().rGetUpperCorner().rGetLocation();
        c_vector<double, 2> lower_corner = p_mesh->CalculateBoundingBox().rGetLowerCorner().rGetLocation();

        c_vector<double, 2> middle_of_mesh = 0.5 * (upper_corner + lower_corner);


#endif // BOOST_VERSION >= 105200
    }

    void TestGenerateCirclularVoronoi() throw(Exception)
    {
    #if BOOST_VERSION >= 105200

        unsigned num_elems_to_keep = 100;
        unsigned approx_rad = std::ceil(sqrt(num_elems_to_keep/M_PI)) + 1;

        // Generate a mesh that is 10 cells wide, 10 high, with 3 Lloyd's relaxation steps and target average element area 0.5
        VoronoiVertexMeshGenerator generator(2 * approx_rad, 2 * approx_rad, 10, 0.5);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMeshAfterReMesh();

        PRINT_VARIABLE(p_mesh->GetNumAllElements());

        c_vector<double, 2> upper_corner = p_mesh->CalculateBoundingBox().rGetUpperCorner().rGetLocation();
        c_vector<double, 2> lower_corner = p_mesh->CalculateBoundingBox().rGetLowerCorner().rGetLocation();

        c_vector<double, 2> middle_of_mesh = 0.5 * (upper_corner + lower_corner);

        std::vector<std::pair<unsigned, double> > elem_to_dist(p_mesh->GetNumAllElements());
        for (unsigned elem_idx = 0; elem_idx < p_mesh->GetNumAllElements(); elem_idx++)
        {
            elem_to_dist[elem_idx].first = elem_idx;
            elem_to_dist[elem_idx].second = norm_2(p_mesh->GetVectorFromAtoB(middle_of_mesh, p_mesh->GetCentroidOfElement(elem_idx)));
        }

        std::sort(elem_to_dist.begin(), elem_to_dist.end(), CustomComparisonForPair);


        for (unsigned i = 0; i < p_mesh->GetNumAllElements() - num_elems_to_keep; i++)
        {
            p_mesh->DeleteElementPriorToReMesh(elem_to_dist[i].first);
        }

        p_mesh->ReMesh();

        // Create a vertex mesh writer
        VertexMeshWriter<2,2> vertex_mesh_writer("TestGenerateCircularVoronoi", "voronoi_circular_50");
        vertex_mesh_writer.WriteFilesUsingMesh(*p_mesh);
        vertex_mesh_writer.WriteVtkUsingMesh(*p_mesh,"");

        TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), num_elems_to_keep);

        // Convert to a potts mesh 

        double lattice_size = 0.1;
        double domain_size = 20;

        double num_nodes_across = domain_size/lattice_size;

        PottsMeshGenerator<2> potts_generator(num_nodes_across, 0, 0, num_nodes_across, 0, 0); 

        PottsMesh<2>* p_potts_mesh = potts_generator.GetMesh();

        p_potts_mesh->Translate(-0.5*num_nodes_across, -0.5*num_nodes_across);
        p_potts_mesh->Scale(lattice_size,lattice_size,lattice_size);    

        GeneratePottsFromVertexMesh(p_potts_mesh,p_mesh);

        PottsMeshWriter<2> potts_mesh_writer("TestGeneratePottsFromVertex", "potts_mesh");

        TS_ASSERT_EQUALS(p_potts_mesh->GetNumAllElements(), num_elems_to_keep);
        TS_ASSERT_EQUALS(p_potts_mesh->GetNumNodes(), num_nodes_across*num_nodes_across);


        // Write and check it's correct
        potts_mesh_writer.WriteFilesUsingMesh(*p_potts_mesh);


    #endif // BOOST_VERSION >= 105200
    }
};

#endif /*TESTGENERATEVORONOIINITIALCONDITIONS_HPP_*/
