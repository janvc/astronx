# Copyright 2012-2015 Jan von Cosel
#
# This file is part of astronx.
#
# astronx if free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# astronx is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have recieved a copy of the GNU General Public License
# along with astronx. If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################
#

FC = gfortran
MKDIR_P = mkdir -p

SRC_DIR = source
OBJ_DIR = objects
MOD_DIR = modules
BIN_DIR = bin

ASTRONX   = $(addprefix $(BIN_DIR)/,astronx)
EXTRACTOR = $(addprefix $(BIN_DIR)/,extractor)
GENERATOR = $(addprefix $(BIN_DIR)/,generator)

CFLAGS = -Wall -Wextra -g --pedantic-errors -std=f2003 -O0 -fall-intrinsics
VPATH = $(SRC_DIR)/astronx:   \
        $(SRC_DIR)/extractor: \
        $(SRC_DIR)/generator: \
        $(SRC_DIR)/common:    \
        $(OBJ_DIR)

ASTRONX_OBJS = $(addprefix $(OBJ_DIR)/,\
                 modules.o        \
                 read_input.o     \
                 common_utils.o   \
                 astronx_utils.o  \
                 bulirsch_stoer.o \
                 propagate.o      \
                 astronx.o        \
                 )
                 
EXTRACT_OBJS = $(addprefix $(OBJ_DIR)/,\
                 modules.o         \
                 common_utils.o    \
                 extractor_utils.o \
                 trjconv.o         \
                 enprop.o          \
                 extractor.o       \
                 )

GENERAT_OBJS = $(addprefix $(OBJ_DIR)/,\
                 common_utils.o \
                 mkdist.o    \
                 generator.o \
                 )

all: directories $(ASTRONX) $(EXTRACTOR) $(GENERATOR)

directories: $(OBJ_DIR) $(MOD_DIR) $(BIN_DIR)

$(GENERATOR): $(GENERAT_OBJS)
	$(FC) -o $@ $^

$(EXTRACTOR): $(EXTRACT_OBJS)
	$(FC) -o $@ $^

$(ASTRONX): $(ASTRONX_OBJS)
	$(FC) -o $@ $^

$(OBJ_DIR)/astronx.o: $(addprefix $(OBJ_DIR)/,\
                        modules.o       \
                        read_input.o    \
                        common_utils.o  \
                        astronx_utils.o \
                        propagate.o     \
                        )

$(OBJ_DIR)/extractor.o: $(addprefix $(OBJ_DIR)/,\
                          modules.o         \
                          extractor_utils.o \
                          trjconv.o         \
                          enprop.o          \
                          )

$(OBJ_DIR)/generator.o: $(addprefix $(OBJ_DIR)/,\
                          modules.o \
                          mkdist.o  \
                          )

$(OBJ_DIR)/propagate.o: $(addprefix $(OBJ_DIR)/,\
                          modules.o        \
                          read_input.o     \
                          astronx_utils.o  \
                          bulirsch_stoer.o \
                          )

$(OBJ_DIR)/bulirsch_stoer.o: $(addprefix $(OBJ_DIR)/,\
                               modules.o       \
                               read_input.o    \
                               common_utils.o  \
                               astronx_utils.o \
                               )

$(OBJ_DIR)/astronx_utils.o: $(addprefix $(OBJ_DIR)/,\
                              modules.o      \
                              read_input.o   \
                              common_utils.o \
                              )

$(OBJ_DIR)/read_input.o: $(addprefix $(OBJ_DIR)/,\
                           modules.o \
                           )

$(OBJ_DIR)/common_utils.o: $(addprefix $(OBJ_DIR)/,\
                             modules.o \
                             )

$(OBJ_DIR)/extractor_utils.o: $(addprefix $(OBJ_DIR)/,\
                                modules.o \
                                )

$(OBJ_DIR)/trjconv.o: $(addprefix $(OBJ_DIR)/,\
                        modules.o         \
                        extractor_utils.o \
                        )

$(OBJ_DIR)/enprop.o: $(addprefix $(OBJ_DIR)/,\
                       modules.o         \
                       extractor_utils.o \
                       )

$(OBJ_DIR)/mkdist.o: $(addprefix $(OBJ_DIR)/,\
                       modules.o      \
                       common_utils.o \
                       )

$(OBJ_DIR)/%.o: %.f03
	$(FC) -c $(CFLAGS) -J$(MOD_DIR) $< -o $@

$(OBJ_DIR):
	$(MKDIR_P) $(OBJ_DIR)
  
$(MOD_DIR):
	$(MKDIR_P) $(MOD_DIR)
  
$(BIN_DIR):
	$(MKDIR_P) $(BIN_DIR)

clean:
	test -d $(BIN_DIR) && rm -r $(OBJ_DIR) $(MOD_DIR) $(BIN_DIR)
