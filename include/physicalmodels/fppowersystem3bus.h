/*! @file fppowersystem3bus.h
    @brief 3-Bus Power System model
    @author Debojyoti Ghosh

  Fokker-Planck model for a 3-Bus power system: This model governs the dynamics of an electrical system with three generators. The third generator is the reference generator and the state of the other two generators are defined with respect to the state of the third generator. \f$p\left({\bf x}\right)\f$ is the probability density function (PDF), governed by the following advection-diffusion equation:
  \f{equation}{
    \frac{\partial p}{\partial t} + \nabla_{\bf x} \cdot \left(\left<{\bf v}\right> p\right) = \nabla_{\bf x} \cdot \left( {\bf D} \nabla_{\bf x} p \right),
  \f}
  where
  \f{align}{
    {\bf x} &= \left[ \theta_1, \theta_2, \omega_1, \omega_2 \right]^T, \\
    \left<{\bf v}\right> &= \left[ \omega_B\omega_1, \omega_B\omega_2, F_1-\gamma\omega_1-S_1, F_2-\gamma\omega_2-S_2  \right]^T \\
    {\bf D} &= \left[ \begin{array}{cc} {\bf 0}_N & {\bf 0}_N \\ {\bf D}^\theta_N & {\bf D}^\omega_N \end{array} \right]
  \f}
  The parameters for the system are:
  \f{align}{
    N&=2,\\
    D_{ij}^\omega &= \sigma_{ij}^2\lambda_{ij}\left(1-\gamma\lambda_{ij}\right), \\
    D_{ij}^\theta &= \sigma_{ij}^2\lambda_{ij}^2\omega_B, \\
    F_i &= \frac{1}{2H_i} \left<P^m_i\right> - \frac{1}{2H_{ref}}\left<P_{ref}^m\right>, \\
    S_i &= \frac{1}{2H_i} P^e_i - \frac{1}{2H_{ref}} P_{ref}^e, \\
    P_i^e &= \mathcal{R}\left\{ E_i \sum_{k=1}^{N+1} E_k^* \left( Y^* \right)_{ik} \right\}, Y = G+iB
  \f}
  where \f$*\f$ denotes the complex conjugate.

  Reference: To be added soon.
*/

/*! Fokker-Planck model for 3-Bus power system */ 
#define _FP_POWER_SYSTEM_3BUS_  "fp-power-system-3bus"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
/*! Number of spatial dimensions */
#define _MODEL_NDIMS_ 4
/*! Number of variables per grid point */
#define _MODEL_NVARS_ 1

/*! \def FPPowerSystem3Bus
    \brief Structure containing variable and parameters specific to the 3-bus power system model.
    This structure contains the physical parameters and variables for the Fokker-Planck model for
    a 3-bus power system.
*/
/*! \brief Structure containing variable and parameters specific to the 3-bus power system model.
    This structure contains the physical parameters and variables for the Fokker-Planck model for
    a 3-bus power system.
*/
typedef struct fp_power_system__3bus_parameters {

  int    N; /*!< Number indepedent buses in power system (excluding reference bus) (must be 2) */
  double  Pm1_avg,          /*!< Average mechanical power of generator 1 \f$\left< P^m_1\right>\f$ */
          Pm2_avg,          /*!< Average mechanical power of generator 2 \f$\left< P^m_2\right>\f$ */
          Pmref_avg,        /*!< Average mechanical power of generator 3 (reference) \f$\left< P^m_{ref}\right>\f$ */
          H1,               /*!< Inertia of generator 1 \f$ H_1 \f$ */
          H2,               /*!< Inertia of generator 2 \f$ H_2 \f$ */
          Href,             /*!< Inertia of generator 3 (reference) \f$ H_{ref} \f$ */
          E1,               /*!< Internal emf phasor for generator 1 \f$ E_1 \f$ */
          E2,               /*!< Internal emf phasor for generator 2 \f$ E_2 \f$ */
          Eref,             /*!< Internal emf phasor for generator 3 (reference) \f$ E_{ref} \f$ */
          omegaB,           /*!< Rotor base angular frequency \f$\omega_B\f$ */
          sigma[2][2],      /*!< Covariance matrix \f$\sigma\f$ */
          lambda[2][2],     /*!< Correlation time \f$\lambda\f$ */
          gamma,            /*!< Damping rate \f$ \gamma \f$ */
          *G,               /*!< Conductance matrix \f$ G \f$ */
          *B;               /*!< Susceptance matrix \f$ B \f$ */

  double pdf_integral;      /*!< Integral of the probability density function over the domain */

} FPPowerSystem3Bus;

int FPPowerSystem3BusInitialize    (void*,void*);
int FPPowerSystem3BusCleanup       (void*);
