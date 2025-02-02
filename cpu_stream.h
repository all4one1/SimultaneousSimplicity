namespace stream_cpu
{
    double dx1(unsigned int l, double* arr) {
        return (arr[l + 1] - arr[l - 1]) / (2.0 * host.hx);
    };

    double dy1(unsigned int l, double* arr) {
        return (arr[l + host.offset] - arr[l - host.offset]) / (2.0 * host.hy);
    };

    double dx2(unsigned int l, double* arr) {
        return (arr[l + 1] - 2.0 * arr[l] + arr[l - 1]) / (host.hx * host.hx);
    };

    double dy2(unsigned int l, double* arr) {
        return (arr[l + host.offset] - 2.0 * arr[l] + arr[l - host.offset]) / (host.hy * host.hy);
    };

    double dx1_eq_0_forward(unsigned int l, double* f) {
        return (4.0 * f[l + 1] - f[l + 2]) / 3.0;
    }
    double dx1_eq_0_back(unsigned int l, double* f) {
        return (4.0 * f[l - 1] - f[l - 2]) / 3.0;
    }
    double dy1_eq_0_up(unsigned int l, double* f) {
        return (4.0 * f[l + host.offset] - f[l + 2 * host.offset]) / 3.0;
    }
    double dy1_eq_0_down(unsigned int l, double* f) {
        return (4.0 * f[l - host.offset] - f[l - 2 * host.offset]) / 3.0;
    }

    #define VX_ dy1(l, ksi)
    #define VY_ -dx1(l, ksi)

    void vorticity(double* omega_new, double* omega, double* ksi, double* T, double* C)
    {
        auto InnerComputing = [&](unsigned int l) {
            return omega[l] + host.tau * (
                (dx1(l, ksi) * dy1(l, omega) - dy1(l, ksi) * dx1(l, omega)) //nonlinear term
                + (dx2(l, omega) + dy2(l, omega)) /** host.Pr*/

                + host.grav_y * host.Ra / host.Pr * (dx1(l, T) - host.density_x)
                - host.grav_x * host.Ra / host.Pr * (dy1(l, T) - host.density_y)

                + host.grav_y * host.Ra / host.Pr * host.K * (dx1(l, C) - host.density_x)
                - host.grav_x * host.Ra / host.Pr * host.K * (dy1(l, C) - host.density_y)
                );
        };

        for (unsigned int j = 0; j <= host.ny; ++j) {
            for (unsigned int i = 0; i <= host.nx; ++i) {
                unsigned int l = i + host.offset * j;

                if (l < host.N) {
                    /* INNER */
                    if (i > 0 && i < host.nx && j > 0 && j < host.ny) {
                        omega_new[l] = InnerComputing(l);
                    }
                    else {
                        if (j == 0 && (i > 0 && i < host.nx)) {
                            omega_new[l] = -0.5 / (host.hy * host.hy) * (8.0 * ksi[l + host.offset] - ksi[l + host.offset * 2]);
                            continue;
                        }
                        else if (j == host.ny && (i > 0 && i < host.nx)) {
                            omega_new[l] = -0.5 / (host.hy * host.hy) * (8.0 * ksi[l - host.offset] - ksi[l - host.offset * 2]);
                            continue;
                        }

                        if (host.xbc == 0) { // closed
                            if (i == 0 && (j > 0 && j < host.ny))
                                omega_new[l] = -0.5 / (host.hx * host.hx) * (8.0 * ksi[l + 1] - ksi[l + 2]);
                            if (i == host.nx && (j > 0 && j < host.ny))
                                omega_new[l] = -0.5 / (host.hx * host.hx) * (8.0 * ksi[l - 1] - ksi[l - 2]);
                            continue;
                        }
                        else if (host.xbc == 1) { // periodic
                            if (i == 0 && (j > 0 && j < host.ny)) {
                                unsigned int ll = host.nx - 1 + host.offset * j;
                                omega_new[l] = InnerComputing(ll);
                                continue;
                            }
                            if (i == host.nx && (j > 0 && j < host.ny)) {
                                unsigned int ll = 1 + host.offset * j;
                                omega_new[l] = InnerComputing(ll);
                                continue;
                            }
                        }
                        {
                            omega_new[l] = 0;
                            omega_new[l] = 0;
                        }
                    }
                }
            }
        }
    }

    void poisson_stream(double* ksi_new, double* ksi, double* omega)
    {
        double tau = 0.2 * host.hx * host.hy; // Вычисляем tau

        for (unsigned int j = 0; j <= host.ny; ++j) {
            for (unsigned int i = 0; i <= host.nx; ++i) {
                unsigned int l = i + host.offset * j;

                if (l < host.N) {
                    /* INNER */
                    if (i > 0 && i < host.nx && j > 0 && j < host.ny) {
                        ksi_new[l] = ksi[l] + tau * (dx2(l, ksi) + dy2(l, ksi) + omega[l]);
                    }
                    else {
                        if (j == 0 && (i > 0 && i < host.nx)) {
                            ksi_new[l] = 0.0;
                            continue;
                        }
                        else if (j == host.ny && (i > 0 && i < host.nx)) {
                            ksi_new[l] = 0.0;
                            continue;
                        }

                        if (host.xbc == 0) { // closed
                            if (i == 0 && (j > 0 && j < host.ny))
                                ksi_new[l] = 0.0;
                            if (i == host.nx && (j > 0 && j < host.ny))
                                ksi_new[l] = 0.0;
                            continue;
                        }
                        else if (host.xbc == 1) { // periodic
                            if (i == 0 && (j > 0 && j < host.ny)) {
                                unsigned int ll = host.nx - 1 + host.offset * j;
                                ksi_new[l] = ksi[ll] + tau * (dx2(ll, ksi) + dy2(ll, ksi));
                                continue;
                            }
                            if (i == host.nx && (j > 0 && j < host.ny)) {
                                unsigned int ll = 1 + host.offset * j;
                                ksi_new[l] = ksi[ll] + tau * (dx2(ll, ksi) + dy2(ll, ksi));
                                continue;
                            }
                        }
                        {
                            ksi_new[l] = 0;
                        }
                    }
                }
            }
        }
    }

    void temperature_2d(double* T, double* T0, double* ksi)
    {
        for (unsigned int j = 0; j <= host.ny; ++j) {
            for (unsigned int i = 0; i <= host.nx; ++i) {
                unsigned int l = i + host.offset * j;

                if (l < host.N) {
                    /* INNER */
                    if (i > 0 && i < host.nx && j > 0 && j < host.ny) {
                        T[l] = T0[l]
                            + host.tau * (
                                -dy1(l, ksi) * dx1(l, T0) + dx1(l, ksi) * dy1(l, T0)
                                + (VX_ * host.density_x + VY_ * host.density_y)
                                + (dx2(l, T0) + dy2(l, T0)) / host.Pr
                                );
                        continue;
                    }
                    else {
                        if (j == 0) {
                            T[l] = 0.0;
                            continue;
                        }
                        else if (j == host.ny) {
                            T[l] = 0.0;
                            continue;
                        }

                        if (host.xbc == 0) { // closed
                            if (i == 0 && (j > 0 && j < host.ny)) {
                                T[l] = dx1_eq_0_forward(l, T0);
                                continue;
                            }
                            if (i == host.nx && (j > 0 && j < host.ny)) {
                                T[l] = dx1_eq_0_back(l, T0);
                                continue;
                            }
                        }
                        else if (host.xbc == 1) { // periodic
                            if (i == 0 && (j > 0 && j < host.ny)) {
                                unsigned int ll = host.nx - 1 + host.offset * j;
                                T[l] = T0[ll];
                                continue;
                            }
                            if (i == host.nx && (j > 0 && j < host.ny)) {
                                unsigned int ll = 1 + host.offset * j;
                                T[l] = T0[ll];
                                continue;
                            }
                        }

                        T[l] = 0;
                    }
                }
            }
        }
    }

    void concentration_2d(double* C, double* C0, double* ksi)
    {
        for (unsigned int j = 0; j <= host.ny; ++j) {
            for (unsigned int i = 0; i <= host.nx; ++i) {
                unsigned int l = i + host.offset * j;

                if (l < host.N) {
                    /* INNER */
                    if (i > 0 && i < host.nx && j > 0 && j < host.ny) {
                        C[l] = C0[l]
                            + host.tau * (
                                -dy1(l, ksi) * dx1(l, C0) + dx1(l, ksi) * dy1(l, C0)
                                + (VX_ * host.density_x + VY_ * host.density_y)
                                + (dx2(l, C0) + dy2(l, C0)) / (host.Le * host.Pr)
                                );
                        continue;
                    }
                    else {
                        if (j == 0) {
                            C[l] = dy1_eq_0_up(l, C0);
                            continue;
                        }
                        else if (j == host.ny) {
                            C[l] = dy1_eq_0_down(l, C0);
                            continue;
                        }

                        if (host.xbc == 0) { // closed
                            if (i == 0 && (j > 0 && j < host.ny)) {
                                C[l] = dx1_eq_0_forward(l, C0);
                                continue;
                            }
                            if (i == host.nx && (j > 0 && j < host.ny)) {
                                C[l] = dx1_eq_0_back(l, C0);
                                continue;
                            }
                        }
                        else if (host.xbc == 1) { // periodic
                            if (i == 0 && (j > 0 && j < host.ny)) {
                                unsigned int ll = host.nx - 1 + host.offset * j;
                                C[l] = C0[ll];
                                continue;
                            }
                            if (i == host.nx && (j > 0 && j < host.ny)) {
                                unsigned int ll = 1 + host.offset * j;
                                C[l] = C0[ll];
                                continue;
                            }
                        }

                        C[l] = 0;
                    }
                }
            }
        }
    }

    void swap_three(double* f_old, double* f_new, double* f2_old, double* f2_new, double* f3_old, double* f3_new)
    {
        for (unsigned int l = 0; l < host.N; l++)
        {
            f_old[l] = f_new[l];
            f2_old[l] = f2_new[l];
            f3_old[l] = f3_new[l];
        }
    }


    struct CuPoisson
    {
        unsigned int k = 0;
        double eps = 0, res = 0, res0 = 0;
        double eps_iter = 1e-4;

        void solve(double *ksi, double *ksi0, double *omega)
        {
            k = 0;
            eps = 1.0;
            res = 0.0;
            res0 = 0.0;

            auto reduce = [this](double* f)
            {
                double s = 0;
                for (unsigned int l = 0; l < host.N; l++)
                    s += abs(f[l]);
                return s;
            };


            for (k = 1; k < 1000000; k++)
            {
                poisson_stream(ksi, ksi0, omega);
                res = reduce(ksi);
                eps = abs(res - res0) / (res0 + 1e-5);
                res0 = res;

                std::swap(ksi, ksi0);

                if (eps < eps_iter)	break;
                if (k % 1000 == 0) std::cout << "device k = " << k << ", eps = " << eps << std::endl;
            }
            if (k > 100) std::cout << "device k = " << k << ", eps = " << eps << std::endl;
        }
    };


}