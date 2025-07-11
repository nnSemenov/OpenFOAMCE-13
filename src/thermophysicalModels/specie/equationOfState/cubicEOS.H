//
// Created by joseph on 25-7-11.
//

#ifndef CUBICEOS_H
#define CUBICEOS_H

#include "autoPtr.H"
#include "mathematicalConstants.H"
#include <cassert>

namespace Foam {

    struct cubicEOSCoefficient {
        scalar a;
        scalar b;
        scalar A;
        scalar B;
        scalar Tr;
        scalar Pr;
        scalar kappa;
        scalar alpha;
    };

    // Solve z^3 + a2 *z^2 + a1 * z + z0 = 0
    scalar solveCubicEquation(scalar a2, scalar a1, scalar a0) {

        const scalar Q = (3*a1 - a2*a2)/9.0;
        const scalar Rl = (9*a2*a1 - 27*a0 - 2*a2*a2*a2)/54.0;

        const scalar Q3 = Q*Q*Q;
        const scalar D = Q3 + Rl*Rl;

        scalar root = -1;

        if (D <= 0)
        {
            const scalar th = ::acos(Rl/sqrt(-Q3));
            const scalar qm = 2*sqrt(-Q);
            const scalar r1 = qm*cos(th/3.0) - a2/3.0;
            const scalar r2 =
                    qm*cos((th + 2*constant::mathematical::pi)/3.0) - a2/3.0;
            const scalar r3 =
                    qm*cos((th + 4*constant::mathematical::pi)/3.0) - a2/3.0;

            root = max(r1, max(r2, r3));
        }
        else
        {
            // One root is real
            const scalar D05 = sqrt(D);
            const scalar S = pow(Rl + D05, 1.0/3.0);
            scalar Tl = 0;
            if (D05 > Rl)
            {
                Tl = -pow(mag(Rl - D05), 1.0/3.0);
            }
            else
            {
                Tl = pow(Rl - D05, 1.0/3.0);
            }

            root = S + Tl - a2/3.0;
        }
        assert(root>0);

        return root;
    }

    struct residualEnergy {
        scalar h_R;//[J/kg]
        scalar e_R;//[J/kg]
        scalar Z;
    };

    struct residualSpecificHeat {
        scalar Cp_R;//[J/kg/K]
        scalar Cv_R;//[J/kg/K]
        scalar Cp_R_minus_Cv_R;//[J/kg/K]
        scalar Z;
    };

}

#endif // CUBICEOS_H
