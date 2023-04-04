"""
Custom ChimeraX commands

See API:
https://www.cgl.ucsf.edu/chimerax/docs/devel/py-modindex.html

guide:
to print to chimerax console, use `session.logger.info("string")

"""


def straighten(
    session,
    residues,
    phi=-165.0,
    psi=165.0,
    omega=155.0,
    noise=3.0,
    move_smaller=True,
):
    import numpy as np
    import random

    nres = len(residues)

    phi_arr = np.random.normal(loc=phi, scale=noise, size=nres)
    psi_arr = np.random.normal(loc=psi, scale=noise, size=nres)

    for i, res in enumerate(residues):
        if res.name == "PRO":
            # res.phi = np.random.normal(loc=-75.0, scale=noise, size=1)[0]
            # res.psi = random.gauss(psi, noise)
            # res.omega = random.gauss(omega, noise)
            res.set_psi(
                random.gauss(psi, noise), move_smaller_side=move_smaller
            )
            res.set_omega(
                random.gauss(omega, noise), move_smaller_side=move_smaller
            )
        else:
            # res.phi = phi_arr[i]
            # res.psi = psi_arr[i]
            res.set_phi(phi_arr[i], move_smaller_side=move_smaller)
            res.set_psi(psi_arr[i], move_smaller_side=move_smaller)


def register_command(session):
    from chimerax.core.commands import CmdDesc, register, FloatArg, BoolArg
    from chimerax.atomic import ResiduesArg

    desc = CmdDesc(
        required=[("residues", ResiduesArg)],
        optional=[
            ("noise", FloatArg),
            ("phi", FloatArg),
            ("psi", FloatArg),
            ("omega", FloatArg),
            ("move_smaller", BoolArg),
        ],
        synopsis="straighten specified residues",
    )

    register("straighten", desc, straighten, logger=session.logger)


register_command(session)
