#########################################################################
# 
#                 #######               ######  #######
#                 ##                    ##   ## ##
#                 #####   ##  ## #####  ##   ## ## ####
#                 ##       ####  ## ##  ##   ## ##   ##
#                 ####### ##  ## ###### ######  #######
#
#  ExaDG - High-Order Discontinuous Galerkin for the Exa-Scale
#
#  Copyright (C) 2021 by the ExaDG authors
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#########################################################################

FUNCTION(SERIALIZE_CMD_COMMAND CMD_COMMAND_STRING_WITH_ARGUMENTS CMD_OUTPUT)

  SEPARATE_ARGUMENTS(CMD_COMMAND UNIX_COMMAND ${CMD_COMMAND_STRING_WITH_ARGUMENTS})
  
  EXECUTE_PROCESS(
    COMMAND ${CMD_COMMAND}
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
    OUTPUT_VARIABLE CMD_OUT
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  SET(${CMD_OUTPUT} ${CMD_OUT} PARENT_SCOPE)

ENDFUNCTION(SERIALIZE_CMD_COMMAND)

MACRO(FETCH_GIT_INFO GIT_BRANCH GIT_SHORTREV)
  SERIALIZE_CMD_COMMAND("git rev-parse --abbrev-ref HEAD" EXADG_GIT_BRANCH)
  SERIALIZE_CMD_COMMAND("git describe --always --tags" EXADG_GIT_SHORTREV)
ENDMACRO(FETCH_GIT_INFO)
