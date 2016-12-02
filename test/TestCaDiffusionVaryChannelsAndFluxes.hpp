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

#ifndef TESTCADIFFUSIONVARYCHANNELSANDFLUXES_HPP_
#define TESTCADIFFUSIONVARYCHANNELSANDFLUXES_HPP_

/*
 * [[PageOutline]]
 *
 * = Ca^2+^ Channel Re-localization to Plasma-Membrane Microdomains Strengthens Activation of Ca^2+^-Dependent Nuclear Gene Expression  =
 *
 * Code to accompany the paper [http://dx.doi.org/10.1016/j.celrep.2015.06.018 Samanta et al. 2015].
 *
 * == Code Walkthrough ==
 *
 * The following wiki page provides a walk-through of the Chaste code
 * that was used to perform the simulations in this paper.
 *
 * First we include some header files:
 */

#include <cxxtest/TestSuite.h>

#include <boost/assign/list_of.hpp>

#include "GmshMeshReader.hpp"
#include "UblasIncludes.hpp"
//#include "SimpleLinearEllipticSolver.hpp"
#include "SimpleLinearParabolicSolver.hpp"
#include "TetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "CommandLineArguments.hpp"

#include "PetscSetupAndFinalize.hpp"

/*
 * === Set up a diffusion equation with a source term ===
 *
 * d[Ca]/dt = D_Ca Laplacian([Ca]) + Q
 *
 * [Ca] in units of uM
 * D_Ca = 300 (nm)^2^/us
 * integral of Q per ion channel's worth of elements over which is to be applied = 2.5133e4 uM / us
 */
template <unsigned SPACE_DIM>
class DiffusionEquationWithSourceTerm : public AbstractLinearParabolicPde<SPACE_DIM>
{
private:
    const std::vector<c_vector<double, SPACE_DIM> >& mrChannelLocations;

    /** The elements that are source elements */
    std::vector<std::vector<unsigned> > mElementsEachChannel;

    std::vector<double> mConcentrationSourcePerUnitVolume; // in units of uM / us

public:
    DiffusionEquationWithSourceTerm(AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>* pMesh,
                                    const std::vector<c_vector<double, SPACE_DIM> >& rChannelLocations)
     : AbstractLinearParabolicPde<SPACE_DIM>(),
       mrChannelLocations(rChannelLocations)
    {
        // We're going to work out the total volume of the elements that are going to have a source term
        // in them here.
        std::vector<double> source_volumes(mrChannelLocations.size(), 0.0);
        mConcentrationSourcePerUnitVolume.resize(mrChannelLocations.size());
        mElementsEachChannel.resize(mrChannelLocations.size());

        // Make a list of elements that we are going to say are source elements.
        for (typename TetrahedralMesh<SPACE_DIM,SPACE_DIM>::ElementIterator elt_iter = pMesh->GetElementIteratorBegin();
             elt_iter != pMesh->GetElementIteratorEnd();
             ++elt_iter)
        {
            // Decide whether this should be classed as a source element or not.
            c_vector<double, SPACE_DIM> location = elt_iter->CalculateCentroid();

            for (unsigned channel=0; channel < mrChannelLocations.size(); channel++)
            {
                if (norm_2(location-mrChannelLocations[channel]) <= 3.0) // If centroid is within 3nm of channel say it is a source.
                {
                    //std::cout << "Source Element Recorded\n";

                    // Make a note of this element
                    mElementsEachChannel[channel].push_back(elt_iter->GetIndex());

                    // And add its volume to the total volume of source.
                    c_matrix<double, 3, 3> jacob;
                    double det;
                    elt_iter->CalculateJacobian(jacob, det);
                    source_volumes[channel] += elt_iter->GetVolume(det);
                    //std::cout << "Total source volume[" << channel << "] = " << source_volumes[channel] << "\n";
                }
            }
        }
        // This number is magical and should not be changed unless ion current changed.
        double total_source_required = 2.5133e4; // Fiddly conversion from ionic current in single channel.
        
        if (CommandLineArguments::Instance()->OptionExists("--low-current"))
		{
			total_source_required /= 1.8; // Another magic number from change in driving voltage from 180mV to 100mV.
		}
        
        for (unsigned channel=0; channel<mrChannelLocations.size(); channel++)
        {
            mConcentrationSourcePerUnitVolume[channel] = total_source_required / source_volumes[channel];
        }
    }

    double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& rPoint,
                             double u,
                             Element<SPACE_DIM,SPACE_DIM>* pElement)
    {
        for (unsigned channel=0; channel<mElementsEachChannel.size(); channel++)
        {
            if (std::find(mElementsEachChannel[channel].begin(), mElementsEachChannel[channel].end(), pElement->GetIndex())
                != mElementsEachChannel[channel].end())
            {
                return mConcentrationSourcePerUnitVolume[channel];
            }
        }
        return 0.0;
    }

    /* The Diffusion constant for calcium is 300 um^2^ / s
     * This is equivalent to
     * 300 (nm)^2^ / us
     */
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rPoint,
                                                                Element<SPACE_DIM,SPACE_DIM>* pElement=NULL)
    {
        const double diffusion_constant = 300; // Units : 3000 (nm)^2/us
        return diffusion_constant*identity_matrix<double>(SPACE_DIM);
    }

    double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& rPoint)
    {
        return 1.0;
    }
};

/*
 * === Test class and method to look at Calcium diffusion ===
 */
class TestCaDiffusion : public CxxTest::TestSuite
{
public:
    /*
     * We will solve
     * du/dt = div(grad u) + u, in 3d, with boundary conditions Ca=0 on the boundary, and initial
     * conditions Ca=0.     *
     *
     * Our units throughout this are:
     * distance : nanometers nm
     * time : microseconds us
     */
    void TestSolvingDiffusionEquationForCalcium() throw(Exception)
    {
    	unsigned scenario = UNSIGNED_UNSET;
		if (!CommandLineArguments::Instance()->OptionExists("--scenario"))
		{
			EXCEPTION("Please enter the command line argument '--scenario' options are 1--6.");
		}
		else
		{
			scenario = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("--scenario");
			if (scenario == 0u || scenario >= 7u)
			{
				EXCEPTION("'--scenario' options should be in the range 1--6.");
			}
		}
		
		/* 
		 *  These scenarios are sketched out here:
		 *
		 *  Channels are arranged like this on a honeycomb pattern:
		 *
		 *  channel 2               channel 1
		 *        --              --
		 *          -- channel 0 --
		 *        --              --
		 *  channel 4               channel 3
		 *
		 *  Scenarios 1 and 2: 4 active channels
		 *   - scenario 1: channel 0 closed; 1-4 open      (occurs 4 times due to symmetry)
		 *   - scenario 2: channel 0, 2-4 open; 1 closed   (occurs once)
		 *
		 *  Scenarios 3-6: 2 active channels
		 *   - scenario 3: 0 and 2 open                    (occurs 4 times due to symmetry)
		 *   - scenario 4: 2 and 4 open                    (occurs twice)
		 *   - scenario 5: 1 and 2 open                    (occurs twice)
		 *   - scenario 6: 2 and 3 open                    (occurs twice)
		 */
		std::map<unsigned, std::list<unsigned> > map_from_scenario_to_included_channels;
		std::list<unsigned> included_channels = boost::assign::list_of(1)(2)(3)(4);
		map_from_scenario_to_included_channels[1] = included_channels;
		included_channels = boost::assign::list_of(0)(2)(3)(4);
		map_from_scenario_to_included_channels[2] = included_channels;
		included_channels = boost::assign::list_of(0)(2);
		map_from_scenario_to_included_channels[3] = included_channels;
		included_channels = boost::assign::list_of(2)(4);
		map_from_scenario_to_included_channels[4] = included_channels;
		included_channels = boost::assign::list_of(1)(2);
		map_from_scenario_to_included_channels[5] = included_channels;
		included_channels = boost::assign::list_of(2)(3);
		map_from_scenario_to_included_channels[6] = included_channels;
		
		TetrahedralMesh<3,3> mesh;
        /*
         * Either a 3D Disc shaped thing, centred at (0,0,0), radius 10, and height 1.5.
         * N.B. We're in slightly odd units of 10 nm or 1e-8 m (!)
         *
         * Note a coarse version of the mesh is provided to test simulations, but the ones in
         * the paper were run on the refined version included here:
         */
        bool disc = true;
        if (disc)
        {
            GmshMeshReader<3,3> gmsh_reader("projects/CaDiffusion/test/meshes/CaDiffusion.msh"); // for accurate solution
            //GmshMeshReader<3,3> gmsh_reader("projects/CaDiffusion/test/meshes/CaDiffusionCoarse.msh"); // for quick estimate
            mesh.ConstructFromMeshReader(gmsh_reader);
        }
        /*
         * Or a square slab of membrane we construct on the fly
         *
         * Create a 20 by 20 by 1.5 mesh in 3D, this time using the {{{ConstructRegularSlabMesh}}}
         * method on the mesh. The first parameter is the cartesian space-step
         * and the other three parameters are the width, height and depth of the mesh.
         */
        else
        {
            mesh.ConstructRegularSlabMesh(0.25, 20.0, 20.0, 1.5);
            mesh.Translate(-10.0,-10.0,0.0); // Centre it in x-y plane at origin
        }
        mesh.Scale(10,10,10); // To get into units of nm, radius 100nm, height 15nm.

        /* Create some ion channel locations */
        double z_location = 15.0; //nm - on the top / outer membrane.
        std::vector<c_vector<double, 3u> > channel_locations;
        
        // Create ion channels on a honeycomb shape.
        {
            double inter_channel_spacing = 15.0; // nm - from MD Cahalan (2015) paper

            // co-ordinates of a pentagon's vertices on a honeycomb grid since hex-symmetry in channel.
            double y = inter_channel_spacing*(cos(M_PI/3.0));
            double x = inter_channel_spacing*(sin(M_PI/3.0));
            
            std::list<unsigned> channels_to_include = map_from_scenario_to_included_channels[scenario];

            c_vector<double, 3u> location;
            location[0] = 0.0;
            location[1] = 0.0;
            location[2] = z_location;
            if (std::find(channels_to_include.begin(), channels_to_include.end(), 0u) != channels_to_include.end())
            {
            	channel_locations.push_back(location);
            }
            

            location[0] = x;
            location[1] = y;
            location[2] = z_location;
            if (std::find(channels_to_include.begin(), channels_to_include.end(), 1u) != channels_to_include.end())
            {
            	channel_locations.push_back(location);
            }

            location[0] = -x;
            location[1] = y;
            location[2] = z_location;
            if (std::find(channels_to_include.begin(), channels_to_include.end(), 2u) != channels_to_include.end())
            {
            	channel_locations.push_back(location);
            }

            location[0] = x;
            location[1] = -y;
            location[2] = z_location;
            if (std::find(channels_to_include.begin(), channels_to_include.end(), 3u) != channels_to_include.end())
            {
            	channel_locations.push_back(location);
            }

            location[0] = -x;
            location[1] = -y;
            location[2] = z_location;
            if (std::find(channels_to_include.begin(), channels_to_include.end(), 4u) != channels_to_include.end())
            {
            	channel_locations.push_back(location);
            }
        }

        /* Create the PDE object (defined above) */
        DiffusionEquationWithSourceTerm<3u> pde(&mesh, channel_locations);

        /* Create a new boundary conditions container and specify u=0.0 on the boundary. */
        BoundaryConditionsContainer<3u,3u,1u> bcc; // Templated over element dim, space dim, problem dim.

        ConstBoundaryCondition<3u>* p_bc_for_Ca = new ConstBoundaryCondition<3u>(0.0);
        for (TetrahedralMesh<3u,3u>::BoundaryNodeIterator node_iter = mesh.GetBoundaryNodeIteratorBegin();
             node_iter != mesh.GetBoundaryNodeIteratorEnd();
             ++node_iter)
        {
            double x = (*node_iter)->rGetLocation()[0];
            double y = (*node_iter)->rGetLocation()[1];
            //double z = (*node_iter)->rGetLocation()[2];

            if (disc)
            {
                // If we are on the rim of the disc
                if (x*x + y*y >= 100*100 - 1e-6)
                {
                    bcc.AddDirichletBoundaryCondition(*node_iter, p_bc_for_Ca, 0);
                }
            }
            else
            {
                // If we are on the edge of a slab
                if (fabs(x+100) < 1e-6 || fabs(x-100) < 1e-6 || fabs(y+100)<1e-6 || fabs(y-100)<1e-6)
                {
                    bcc.AddDirichletBoundaryCondition(*node_iter, p_bc_for_Ca, 0);
                }
            }
        }

        SimpleLinearParabolicSolver<3,3> solver(&mesh, &pde, &bcc);

        /* For parabolic problems, initial conditions are also needed. The solver will expect
         * a PETSc vector, where the i-th entry is the initial solution at node i, to be passed
         * in. To create this PETSc {{{Vec}}}, we will use a helper function in the {{{PetscTools}}}
         * class to create a {{{Vec}}} of size num_nodes, with each entry set to 0.0. Then we
         * set the initial condition on the solver. */
        Vec initial_condition = PetscTools::CreateAndSetVec(mesh.GetNumNodes(), 0.0);
        solver.SetInitialCondition(initial_condition);

        /* Next define the start time, end time, and timestep, and set them. */
        double t_start = 0; // micro seconds
        double t_end = 2; // micro seconds
        
        // For production runs
        //double dt = 0.001; // micro seconds
        
        double dt = 0.01; // micro seconds
        solver.SetTimes(t_start, t_end);
        solver.SetTimeStep(dt);

		std::stringstream output_folder;
		output_folder << "CaDiffusion/time_dependent_scenario_" << scenario << "_current_";
		
		if (CommandLineArguments::Instance()->OptionExists("--low-current"))
		{
			output_folder << "low";
		}
		else
		{
			output_folder << "normal";
		}
		
        solver.SetOutputDirectoryAndPrefix(output_folder.str(),"results");
        solver.SetOutputToVtk(true);
        Vec result = solver.Solve();

        // Write a copy of the mesh to examine in a different format.
        TrianglesMeshWriter<3,3> writer("CaDiffusion/mesh", "disc", false);
        writer.WriteFilesUsingMesh(mesh);

        /* All PETSc vectors should be destroyed when they are no longer needed. */
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);        
    }
};

#endif // TESTCADIFFUSIONVARYCHANNELSANDFLUXES_HPP_
