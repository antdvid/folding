#include "bending.h"

double* BendingForce::getExternalForce(SpringVertex* sv) 
{
    	return ((POINT*)(sv->org_vtx))->force;
}

void BendingForce::computeExternalForce()
{
        SURFACE **surf;
        TRI *tri;

        intfc_surface_loop(intfc, surf)
        {
            if (wave_type(*surf) != ELASTIC_BOUNDARY) continue;
            if (is_bdry(*surf)) continue;
            clear_surf_point_force(*surf);
            surf_tri_loop(*surf, tri)
            {
                for (int i = 0; i < 3; ++i)
                {
                    POINT *p1 = Point_of_tri(tri)[i];
                    TRI* n_tri = Tri_on_side(tri, (i+1)%3);
                    if (n_tri != NULL)
                        calculateBendingForce3d2006(p1, tri, n_tri);
                        //calculateBendingForce3d2003(p1, tri, n_tri);
                }
            }
        }
}       /* setBendingForce3d */

void BendingForce::calculateBendingForce3d2003(
        POINT* p1,
        TRI* t1,
        TRI* t2)
{
        if (Mag3d(Tri_normal(t1)) < MACH_EPS ||
            Mag3d(Tri_normal(t2)) < MACH_EPS) return;

        int index = Vertex_of_point(t1, p1);
        POINT *p3 = Point_of_tri(t1)[(index+1)%3];
        POINT *p4 = Point_of_tri(t1)[(index+2)%3];
        index = 3 - Vertex_of_point(t2, p3) - Vertex_of_point(t2, p4);
        POINT *p2 = Point_of_tri(t2)[index];

        double x13[3], x14[3], x23[3], x24[3], E[3];
        for (int i = 0; i < 3; ++i)
        {
            x13[i] = Coords(p1)[i] - Coords(p3)[i];
            x14[i] = Coords(p1)[i] - Coords(p4)[i];
            x23[i] = Coords(p2)[i] - Coords(p3)[i];
            x24[i] = Coords(p2)[i] - Coords(p4)[i];
            E[i] = Coords(p4)[i] - Coords(p3)[i];
        }
        double N1[3], N2[3], E_mag, N1_mag, N2_mag;
        double n1[3], n2[3];
        Cross3d(x13, x14, N1);
        Cross3d(x24, x23, N2);
        E_mag = Mag3d(E);
        N1_mag = Mag3d(N1);
        N2_mag = Mag3d(N2);
        for (int i = 0; i < 3; ++i)
        {
            E[i] = E[i] / E_mag;
            n1[i] = N1[i] / N1_mag;
            n2[i] = N2[i] / N2_mag;
            N1[i] = n1[i] / N1_mag;
            N2[i] = n2[i] / N2_mag;
        }
        double u1[3], u2[3], u3[3], u4[3];
        for (int i = 0; i < 3; ++i)
        {
            u1[i] = E_mag * N1[i];
            u2[i] = E_mag * N2[i];
            u3[i] = Dot3d(x14, E) * N1[i] + Dot3d(x24, E) * N2[i];
            u4[i] = -Dot3d(x13, E) * N1[i] - Dot3d(x23, E) * N2[i];
        }
        double bend_stiff = 0.00001;
        double coeff = bend_stiff * sqr(E_mag) / (N1_mag + N2_mag);
        if (Dot3d(n1, n2) > 1.0 + 1.0e-10)
        {
            printf("t1 = %20.14f %20.14f %20.14f\n", Tri_normal(t1)[0],
                        Tri_normal(t1)[1],Tri_normal(t1)[2]);
            printf("n1 = %20.14f %20.14f %20.14f\n", n1[0],n1[1],n1[2]);
            printf("t2 = %20.14f %20.14f %20.14f\n", Tri_normal(t2)[0],
                        Tri_normal(t2)[1],Tri_normal(t2)[2]);
            printf("n2 = %20.14f %20.14f %20.14f\n", n2[0],n2[1],n2[2]);
            printf("Dot3d(n1, n2) = %20.14f \n", Dot3d(n1, n2));
            printf("u1 = %20.14f %20.14f %20.14f\n", u1[0], u1[1], u1[2]);
            printf("u2 = %20.14f %20.14f %20.14f\n", u2[0], u2[1], u2[2]);
            printf("u3 = %20.14f %20.14f %20.14f\n", u3[0], u3[1], u3[2]);
            printf("u4 = %20.14f %20.14f %20.14f\n", u4[0], u4[1], u4[2]);
            printf("N1 = %20.14f %20.14f %20.14f\n", N1[0], N1[1], N1[2]);
            printf("N2 = %20.14f %20.14f %20.14f\n", N2[0], N2[1], N2[2]);
            clean_up(0);
        }
        double sine_half_theta = sqrt(0.5 * std::max(0.0, 1.0 - Dot3d(n1, n2)));
        double tmp[3];
        Cross3d(n1, n2, tmp);
        if (Dot3d(tmp, E) < 0)
            sine_half_theta *= -1.0;
        coeff *= sine_half_theta;
        double bend_damp = 0.0;
        double dtheta = 0.0;
        dtheta = Dot3d(u1, p1->vel) + Dot3d(u2, p1->vel) +
                 Dot3d(u3, p3->vel) + Dot3d(u4, p4->vel);
        if (fabs(dtheta) < 1.0e-10) dtheta = 0.0;
        if (dtheta > 0.0)
        {
            double sum[3];
            for (int i = 0; i < 3; ++i)
                sum[i] = u1[i] + u2[i] + u3[i] + u4[i];
            printf("sum = %20.14f %20.14f %20.14f\n",sum[0],sum[1],sum[2]);
            printf("p1->vel = %f %f %f\n",p1->vel[0],p1->vel[1],p1->vel[2]);
            printf("dtheta = %20.14f\n", dtheta);
            clean_up(0);
        }
        coeff += -bend_damp * E_mag * dtheta;
        for (int i = 0; i < 3; ++i)
        {
            // each tri_pair will be calculated twice
            p1->force[i] += coeff * u1[i] * 0.5;
            p2->force[i] += coeff * u2[i] * 0.5;
            p3->force[i] += coeff * u3[i] * 0.5;
            p4->force[i] += coeff * u4[i] * 0.5;
        }
}       /* calculateBendingForce3d */

void BendingForce::calculateBendingForce3d2006(
        POINT* p1,
        TRI* t1,
        TRI* t2)
{
        if (Mag3d(Tri_normal(t1)) < MACH_EPS ||
            Mag3d(Tri_normal(t2)) < MACH_EPS) return;

        int index = Vertex_of_point(t1, p1);
        POINT *p3 = Point_of_tri(t1)[(index+1)%3];
        POINT *p4 = Point_of_tri(t1)[(index+2)%3];
        index = 3 - Vertex_of_point(t2, p3) - Vertex_of_point(t2, p4);
        POINT *p2 = Point_of_tri(t2)[index];

        double x13[3], x14[3], x23[3], x24[3], E[3];
        for (int i = 0; i < 3; ++i)
        {
            x13[i] = Coords(p1)[i] - Coords(p3)[i];
            x14[i] = Coords(p1)[i] - Coords(p4)[i];
            x23[i] = Coords(p2)[i] - Coords(p3)[i];
            x24[i] = Coords(p2)[i] - Coords(p4)[i];
            E[i] = Coords(p4)[i] - Coords(p3)[i];
        }
        double N1[3], N2[3], N3[3], N4[3];
        Cross3d(x13, x14, N1);
        Cross3d(x24, x23, N2);
        Cross3d(x23, x13, N3);
        Cross3d(x14, x24, N4);
        double N1_mag, N2_mag, N3_mag, N4_mag;
        N1_mag = Mag3d(N1);
        N2_mag = Mag3d(N2);
        N3_mag = Mag3d(N3);
        N4_mag = Mag3d(N4);
        double a1, a2, a3, a4;
        a1 = N2_mag / (N1_mag + N2_mag);
        a2 = 1 - a1;
        a3 = - N4_mag / (N3_mag + N4_mag);
        a4 = -1 - a3;
        double R[3];
        for (int i = 0; i < 3; ++i)
            R[i] = a1 * Coords(p1)[i] + a2 * Coords(p2)[i] +
                   a3 * Coords(p3)[i] + a4 * Coords(p4)[i];
        double E_mag, h1, h2;
        E_mag = Mag3d(E);
        h1 = fabs(Dot3d(x13, E) / E_mag);
        h1 = sqrt(Dot3d(x13, x13) - h1 * h1);
        h2 = fabs(Dot3d(x23, E) / E_mag);
        h2 = sqrt(Dot3d(x23, x23) - h2 * h2);
        double bend_stiff = 0.001;
        double lambda = 2.0*(h1 + h2)*E_mag*bend_stiff / (3.0*h1*h2*h1*h2);
        for (int i = 0; i < 3; ++i)
        {
            p1->force[i] += -lambda * a1 * R[i] * 0.5;
            p2->force[i] += -lambda * a2 * R[i] * 0.5;
            p3->force[i] += -lambda * a3 * R[i] * 0.5;
            p4->force[i] += -lambda * a4 * R[i] * 0.5;
        }
        if (0) //Mag3d(p1->force) == 0)
        {
        printf("f1 = [%f %f %f]\n", p1->force[0], p1->force[1], p1->force[2]);
        printf("f2 = [%f %f %f]\n", p2->force[0], p2->force[1], p2->force[2]);
        printf("f3 = [%f %f %f]\n", p3->force[0], p3->force[1], p3->force[2]);
        printf("f4 = [%f %f %f]\n", p4->force[0], p4->force[1], p4->force[2]);
        printf("a = [%f %f %f %f]\n\n", a1, a2, a3, a4);
        printf("h1 = %e, h2 = %e\n", h1, h2);
        }
}       /* end calculateBendingForce3d */

void BendingForce::clear_surf_point_force(SURFACE* surf)
{
        TRI* tri;
        surf_tri_loop(surf, tri)
        {
            for (int i = 0; i < 3; ++i)
            {
                POINT *p = Point_of_tri(tri)[i];
                for (int j = 0; j < 3; ++j)
                    p->force[i] = 0.0;
            }
        }
}       /* clear_surf_point_force */

