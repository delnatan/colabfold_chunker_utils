"""
usage:
# report selected residues
sel #1/A:20-50

reportresidues sel

"""


def reportresidues(session, residues):
    from chimerax.core.commands import run

    errors = []

    for i, res in enumerate(residues):
        if res.polymer_type == res.PT_AMINO:
            structure = res.structure
            ch_id = res.chain_id
            resaa = res.get_one_letter_code()
            resnum = res.number
            aa_str = "#%d/%s:%d (%s)" % (structure.id[0], ch_id, resnum, resaa)
            session.logger.info(aa_str)
        else:
            session.logger.info(str(res))

    # form chimera selection syntax

    if errors:
        msg = "\n".join(errors)
        session.logger.warning(msg)


def register_command(session):
    from chimerax.atomic import ResiduesArg
    from chimerax.core.commands import CmdDesc, register

    desc = CmdDesc(
        required=[("residues", ResiduesArg)],
        synopsis="Print selected residues to console",
    )

    register("reportresidues", desc, reportresidues, logger=session.logger)


register_command(session)
