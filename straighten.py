def straighten(session, mobile_sel, mobile_terminal, noise=5.0):
    """
    Sets phi (ϕ) and psi (ψ) angles of specified residues

      ϕ     ψ     ω
    N----Ca----C-----N
         |
         Cb

    Notes
    -----
    When selection is closer to C-terminus
    N<->Ca: N terminal side is longer, smaller_side == "CA"
    Ca<->C: C terminal side is shorter, smaller_side == "C"

    When selection is closer to N-terminus
    N<->Ca: N terminal side is shorter, smaller_side == "N"
    Ca<->C: C terminal side is longer, smaller_side == "CA"

    Args:
    mobile_sel (selection): selection of residue regions to be moved, must be a
    single contiguous region
    mobile_terminal (str): 'N' or 'C' terminal side
    noise (float): gaussian standard deviation to 'jitter' angles
    """
    import random

    mobile_terminal = mobile_terminal.upper()

    assert mobile_terminal.lower() in (
        "c",
        "n",
    ), "mobile_terminal can only be N or C"

    for i, res in enumerate(mobile_sel):
        # find backbone C-alpha
        atom = res.find_atom("CA")

        # list containing 'small side' atom w.r.t C-alpha
        small_side_atoms = []

        # handle proline differently (for now)
        if res.name != "PRO":
            for nb, bond in zip(atom.neighbors, atom.bonds):
                if nb.name == "N":
                    small_side_atoms.append(bond.smaller_side.name)
                elif nb.name == "C":
                    small_side_atoms.append(bond.smaller_side.name)
            # if we want to move the 'N' terminus
            if mobile_terminal == "N":
                # residue is near N-term (phi & psi is the smaller side)
                move_smaller = set(small_side_atoms) == {"N", "CA"}
            elif mobile_terminal == "C":
                # residue is near C-term (phi & psi is the smaller side)
                move_smaller = set(small_side_atoms) == {"CA", "C"}

            psi_assgn = random.normalvariate(mu=135.0, sigma=noise)
            res.set_psi(psi_assgn, move_smaller_side=move_smaller)

            phi_assgn = random.normalvariate(mu=-139.0, sigma=noise)
            res.set_phi(phi_assgn, move_smaller_side=move_smaller)


def register_command(session):
    from chimerax.atomic import ResiduesArg
    from chimerax.core.commands import (
        CmdDesc,
        FloatArg,
        StringArg,
        register,
    )

    desc = CmdDesc(
        required=[("mobile_sel", ResiduesArg), ("mobile_terminal", StringArg)],
        optional=[
            ("noise", FloatArg),
        ],
        synopsis=(
            "straighten specified residues. Choose N/C terminal"
            " side to be moved"
        ),
    )

    register("straighten", desc, straighten, logger=session.logger)


register_command(session)
