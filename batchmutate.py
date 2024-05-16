"""
usage:
# select all disordered residues
select #1 & coil

batchmutate sel
"""


def batchmutate(session, residues):
    from chimerax.core.commands import run

    errors = []
    for i, res in enumerate(residues):
        if res.name == "PRO":
            command_str = "swapaa %s %s log false" % (
                res.string(style="command"),
                "ala",
            )
            try:
                run(session, command_str)
            except (UserError, LimitationError) as e:
                errors.append(str(e))
            session.logger.status(f"{res} -> Ala")

    if errors:
        msg = "\n".join(errors)
        session.logger.warning(msg)


def register_command(session):
    from chimerax.core.commands import CmdDesc, register, FloatArg, BoolArg
    from chimerax.atomic import ResiduesArg

    desc = CmdDesc(
        required=[("residues", ResiduesArg)],
        synopsis="batch mutate prolines to alanines",
    )

    register("batchmutate", desc, batchmutate, logger=session.logger)


register_command(session)
