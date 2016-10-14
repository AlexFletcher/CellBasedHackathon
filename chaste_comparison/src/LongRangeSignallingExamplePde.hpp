
#ifndef LONGRANGESIGNALLINGEXAMPLEPDE_HPP_
#define LONGRANGESIGNALLINGEXAMPLEPDE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellPopulation.hpp"
#include "AbstractLinearParabolicPde.hpp"

/**
 * A parabolic PDE to be solved numerically using the finite element method, for
 * coupling to the long-range signalling example.
 *
 * The PDE takes the form
 *
 * du/dt = Grad.(D*Grad(u)) + k*u*rho(x),
 *
 * where the scalars D and k are specified by the members mDiffusionCoefficient and
 * mSourceCoefficient, respectively. Their values must be set in the constructor.
 *
 * For a node of the finite element mesh with location x, the function rho(x)
 * equals one if there is a non-apoptotic cell associated with x, and
 * zero otherwise. Here, 'associated with' takes a different meaning for each
 * cell population class, and is encoded in the method IsPdeNodeAssociatedWithNonApoptoticCell().
 */
template<unsigned DIM>
class LongRangeSignallingExamplePde : public AbstractLinearParabolicPde<DIM,DIM>
{
private:

    /** Needed for serialization.*/
    friend class boost::serialization::access;
    /**
     * Serialize the PDE and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
       archive & boost::serialization::base_object<AbstractLinearParabolicPde<DIM, DIM> >(*this);
       archive & mDuDtCoefficient;
       archive & mDiffusionCoefficient;
       archive & mSourceCoefficient;
       archive & mSinkCoefficient;
    }

protected:

    /** The cell population member. */
    AbstractCellPopulation<DIM, DIM>& mrCellPopulation;

    /** Coefficient of rate of change term.  */
    double mDuDtCoefficient;

    /** Diffusion coefficient. */
    double mDiffusionCoefficient;

    /** Coefficient of the rate of release within labelled cells. */
    double mSourceCoefficient;

    /** Coefficient of the rate of uptake of the dependent variable within any cell.
     * Note: This variable is to be defined positive, the minus sign is added in code.
     */
    double mSinkCoefficient;

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation reference to the cell population
     * @param duDtCoefficient rate of reaction (defaults to 1.0)
     * @param diffusionCoefficient rate of diffusion (defaults to 1.0)
     * @param sourceCoefficient the source term coefficient (defaults to 0.0)
     */
    LongRangeSignallingExamplePde(AbstractCellPopulation<DIM, DIM>& rCellPopulation,
                               double duDtCoefficient=1.0,
                               double diffusionCoefficient=1.0,
                               double sourceCoefficient=0.0,
                               double sinkCoefficient=0.0);

    /**
     * @return const reference to the cell population (used in archiving).
     */
    const AbstractCellPopulation<DIM>& rGetCellPopulation() const;

    /**
     * Overridden ComputeDuDtCoefficientFunction() method.
     *
     * @return the function c(x) in "c(x) du/dt = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)"
     *
     * @param rX the point in space at which the function c is computed
     */
    virtual double ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& rX);

    /**
     * Overridden ComputeSourceTerm() method. That is never called.
     *
     * @return computed source term.
     *
     * @param rX the point in space at which the nonlinear source term is computed
     * @param u the value of the dependent variable at the point
     * @param pElement the mesh element that x is contained in (optional; defaults to NULL).
     */
    virtual double ComputeSourceTerm(const ChastePoint<DIM>& rX,
                                     double u,
                                     Element<DIM,DIM>* pElement=NULL);

    /**
     * Overridden ComputeSourceTermAtNode() method.
     *
     * Note that for CellWise Parabolic PDEs used with CellBasedParabolicPdeSolver
     * this method returns the coefficient of the linear component of the source term.
     *
     * @return computed source term at a node.
     *
     * @param rNode the node at which the nonlinear source term is computed
     * @param u the value of the dependent variable at the node
     */
    virtual double ComputeSourceTermAtNode(const Node<DIM>& rNode, double u);

    /**
     * Overridden ComputeDiffusionTerm() method.
     *
     * @param rX the point in space at which the diffusion term is computed
     * @param pElement the mesh element that x is contained in (optional; defaults to NULL).
     *
     * @return a matrix.
     */
    virtual c_matrix<double,DIM,DIM> ComputeDiffusionTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement=NULL);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LongRangeSignallingExamplePde)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a LongRangeSignallingExamplePde.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const LongRangeSignallingExamplePde<DIM>* t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM, DIM>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a LongRangeSignallingExamplePde.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, LongRangeSignallingExamplePde<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM, DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)LongRangeSignallingExamplePde<DIM>(*p_cell_population);
}
}
} // namespace ...

#endif /*LONGRANGESIGNALLINGEXAMPLEPDE_HPP_*/
