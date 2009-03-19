#ifndef PHYSICAL_CONSTANTS_H
#define PHYSICAL_CONSTANTS_H

#include <math.h>

#define SI_UNITS 0
#define PLANCK_UNITS 1

#if SI_UNITS

#define GRAVITATIONAL_CONSTANT (1.0F)

#define ELECTRON_MASS          (510.9989F)      /* keV c^-2 */
#define PROTON_MASS            (938.2722F)      /* MeV c^-2 */
#define NEUTRON_MASS           (938.2722F)      /* MeV c^-2 */

#define PLANKS_CONSTANT_H      (4.135667e-15)   /* eV s */
#define PLANKS_CONSTANT_HBAR   (6.582119e-16)   /* eV s */

#define SPEED_OF_LIGHT         (1.0F)

#define PERMITTIVITY_OF_VACUUM (8.854e-12)      /* C^2 N^-1 m^-2 */

#define FINE_STRUCTURE_CONSTANT (1/137.)
#define COULOMB_CONSTANT        (1/(4 * M_PI * PERMITTIVITY_OF_VACUUM))

#endif

#if PLANCK_UNITS

#define SPEED_OF_LIGHT              (1.0)
#define GRAVITATIONAL_CONSTANT      (1.0)
#define PLANKS_CONSTANT_HBAR        (1.0)
#define COULOMB_CONSTANT            (1.0)
#define FINE_STRUCTURE_CONSTANT     (1/137.03599911)

#define ELECTRON_CHARGE             (-1)
#define PROTON_CHARGE               (+1)
#define NEUTRON_CHARGE              (0)

#define PLANKS_CONSTANT_H           ((2 * M_PI) * PLANKS_CONSTANT_HBAR)
#define PERMITTIVITY_OF_VACUUM      (1 / (4 * M_PI))
#define BOLTZMANN_CONSTANT          (0)

#define ELECTRON_MASS               (4.1854e-23)
#define PROTON_MASS                 (7.6849e-20)
#define NEUTRON_MASS                (7.6956e-20)

#define SOLAR_MASS                  (9.13873e37)

#define SECONDS(s)                  (.18548e44 * (s))

#endif

#if !(SI_UNITS || PLANCK_UNITS)
#error "No units defined!"
#endif

#endif
