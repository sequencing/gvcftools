#
# IDS must be defined:
#
IDS := CXXFLAGS LIBUTIL_DIR

assert-defined-indirect = $(if $($1),,$(error Variable '$1' must be defined))
$(foreach I,$(IDS), $(call assert-defined-indirect,$I))

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:%.cpp=%.o)

LIBLABEL=$(notdir $(CURDIR))
LIBNAME=$(LIBLABEL).a

$(OBJS) : CXXFLAGS += $(LIBUTIL_DIR)

$(LIBNAME): $(OBJS)
	$(AR) -csru $@ $(OBJS)

test:

clean:
	rm -f $(LIBNAME) $(OBJS)

