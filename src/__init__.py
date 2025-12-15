# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.toolshed import BundleAPI

# Subclass from chimerax.core.toolshed.BundleAPI and
# override the method for registering commands,
# inheriting all other methods from the base class.
class _MyAPI(BundleAPI):

    api_version = 1     # start_tool called with BundleInfo and
                        # ToolInfo instance (vs. BundleInfo and
                        # tool name when api_version==0 [the default])

    # Override method
    @staticmethod
    def start_tool(session, bi, ti):
        # session is an instance of chimerax.core.session.Session
        # bi is an instance of chimerax.core.toolshed.BundleInfo
        # ti is an instance of chimerax.core.toolshed.ToolInfo

        # This method is called once for each time the tool is invoked.

        # We check the name of the tool, which should match one of the
        # ones listed in bundle_info.xml (without the leading and
        # trailing whitespace), and create and return an instance of the
        # appropriate class from the ``tool`` module.
        if ti.name == "Cilia Builder":
            from . import tool
            return tool.CiliaBuilder(session, ti.name)
        raise ValueError("trying to start unknown tool: %s" % ti.name)

    @staticmethod
    def register_command(bi, ci, logger):
        # bi is an instance of chimerax.core.toolshed.BundleInfo
        # ci is an instance of chimerax.core.toolshed.CommandInfo
        # logger is an instance of chimerax.core.logger.Logger

        # This method is called once for each command listed
        # in bundle_info.xml.  Since we only listed one command,
        # we expect only a single call to this method.

        # We import the function to call and its argument
        # description from the ``cmd`` module, adding a
        # synopsis from bundle_info.xml if none is supplied
        # by the code.
        from . import cmd
        from chimerax.core.commands import register
        
        from . import cmd
        from chimerax.core.commands import register
        
        if ci.name == "ciliabuild":
            func = cmd.ciliabuild
            desc = cmd.ciliabuild_desc
        elif ci.name == "centriolebuild":
            func = cmd.centriolebuild
            desc = cmd.centriolebuild_desc
        else:
            logger.warning(f"Unknown command {ci.name} listed in bundle_info.xml")
            return

        if desc.synopsis is None:
            desc.synopsis = ci.synopsis
            
        # We then register the function as the command callback
        # with the chimerax.core.commands module.
        register(ci.name, desc, func)  
          
# Create the ``bundle_api`` object that ChimeraX expects.
bundle_api = _MyAPI()