"""
Custom ChimeraX commands

See API:
https://www.cgl.ucsf.edu/chimerax/docs/devel/py-modindex.html

guide:
to print to chimerax console, use `session.logger.info("string")

"""


def straighten(
    session, model_id, start_residue, end_residue, chain_id="A", noise=3.0
):
    import numpy as np

    session.logger.info(
        "Straightening #%d: %d to %d" % (model_id, start_residue, end_residue)
    )

    model_id_str = "%d" % (model_id)

    model = [m for m in session.models if m.id_string == model_id_str][0]
    chain = [c for c in model.chains if c.chain_id == chain_id][0]
    # missing residues are 'None', so we dont want them
    residues = [
        chain.residues[i]
        for i in range(start_residue - 1, end_residue)
        if chain.residues[i] is not None
    ]
    nres = len(residues)
    phi_arr = np.random.normal(loc=-140.0, scale=noise, size=nres)
    psi_arr = np.random.normal(loc=130.0, scale=noise, size=nres)
    for i, res in enumerate(residues):
        if res.name == "PRO":
            # res.phi = np.random.normal(loc=-75.0, scale=noise, size=1)[0]
            res.psi = np.random.normal(loc=150, scale=noise, size=1)[0]
            res.omega = np.random.normal(loc=155, scale=noise, size=1)[0]
        else:
            res.phi = phi_arr[i]
            res.psi = psi_arr[i]
        # session.logger.info(
        #     "%s%d; phi=%.1f, psi=%.1f"
        #     % (res.name, res.number, res.phi, res.psi)
        # )


def register_command(session):
    from chimerax.core.commands import (
        CmdDesc,
        register,
        StringArg,
        IntArg,
        FloatArg,
        OpenFolderNameArg,
        BoolArg,
    )

    desc = CmdDesc(
        required=[
            ("model_id", IntArg),
            ("start_residue", IntArg),
            ("end_residue", IntArg),
        ],
        optional=[("chain_id", StringArg), ("noise", FloatArg)],
        synopsis="straighten specified residues",
    )

    register("straighten", desc, straighten, logger=session.logger)


register_command(session)
