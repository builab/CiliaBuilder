# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.commands import CmdDesc
from chimerax.core.commands import FloatArg
from chimerax.core.commands import StringArg

from .draw import draw_doublet
from .draw import draw_cp

def ciliasim(session, length=3500):
    # All command functions are invoked with ``session`` as its
    # first argument.  Useful session attributes include:
    #   logger: chimerax.core.logger.Logger instance
    #   models: chimerax.core.models.Models instance
    draw_doublet(session, length=length)
    draw_cp(session, length=length)

    session.logger.info("Cilia model generated!")

# CmdDesc contains the command description.  For the
# "hello" command, we expect no arguments.
ciliasim_desc = CmdDesc(keyword=[("length", FloatArg)], synopsis = 'Draw microtubule doublet (A and B tubules)')