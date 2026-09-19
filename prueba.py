def simular_imagen(N, inclinacion_grados):
    # Paso 1: generar posiciones iniciales de los rayos
    theta_0_sim = np.random.uniform(0, 2 * np.pi, N)
    r_0_sim = R * np.sqrt(np.random.uniform(0, 1, N))
    y0_sim = r_0_sim * np.cos(theta_0_sim)
    z0_sim = r_0_sim * np.sin(theta_0_sim)

    # Paso 2: integrar cada rayo en el plano xy
    trayectorias_sim = []

    for theta, r in zip(theta_0_sim, r_0_sim):
        ystart = r
        [n_ini, nx_ini, ny_ini] = refindex(xstart, ystart)
        initial_conditions_for_ode = [xstart, ystart, n_ini, 0.0]

        sol = odeint(flow_deriv, initial_conditions_for_ode, tspan)
        x_rayo = sol[:, 0]
        y_rayo = sol[:, 1]

        x_cortado = []
        y_cortado = []
        cae_bh = False

        for i in range(len(x_rayo)):
            px = x_rayo[i]
            py = y_rayo[i]

            if np.sqrt(px**2 + py**2) < A:
                cae_bh = True
                break

            if abs(px) > xmax or abs(py) > ymax:
                break

            x_cortado.append(px)
            y_cortado.append(py)

        trayectorias_sim.append((theta, x_cortado, y_cortado, cae_bh))

    # Paso 3: inclinar el disco al angulo pedido (reutiliza puntos_disco de la parte 3)
    rotacion_disco = R_rot.from_euler('y', inclinacion_grados, degrees=True)
    disco_inclinado_sim = rotacion_disco.apply(puntos_disco)

    # Construimos el arbol de busqueda una sola vez para este disco inclinado
    arbol_disco = cKDTree(disco_inclinado_sim)

    # Paso 4: rotar cada rayo y verificar interseccion con el disco
    y_hits_sim, z_hits_sim, radios_hits_sim = [], [], []

    for i in range(len(trayectorias_sim)):
        theta, x_cortado, y_cortado, cae_bh = trayectorias_sim[i]

        if len(x_cortado) == 0:
            continue

        Rot = Rx(theta)
        puntos_3d = []
        for x_p, y_p in zip(x_cortado, y_cortado):
            punto_plano = np.array([x_p, y_p, 0])
            puntos_3d.append(Rot @ punto_plano)
        puntos_3d = np.array(puntos_3d)

        # Distancia al punto mas cercano del disco, para cada punto del rayo
        distancias_min, _ = arbol_disco.query(puntos_3d)

        for j in range(len(distancias_min)):
            if distancias_min[j] < d_ave:
                px, py, pz = puntos_3d[j]
                radio_interseccion = np.sqrt(px**2 + py**2 + pz**2)
                y_hits_sim.append(y0_sim[i])
                z_hits_sim.append(z0_sim[i])
                radios_hits_sim.append(radio_interseccion)
                break

    return y_hits_sim, z_hits_sim, radios_hits_sim


def graficar_imagen_bh(y_hits, z_hits, radios_hits, titulo):
    plt.figure(figsize=(7, 7))
    plt.gca().set_facecolor('black')
    plt.gcf().set_facecolor('black')

    plt.scatter(y_hits, z_hits, c=radios_hits, cmap='magma', s=4, vmin=Rin, vmax=Rout)

    plt.title(titulo, color='white')
    plt.xlim(-R-5, R+5)
    plt.ylim(-R-5, R+5)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis('off')
    plt.show()