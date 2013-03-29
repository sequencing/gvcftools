#
# IDS must be defined:
#
IDS := CXXFLAGS LIBTRIO LIBUTIL_DIR

assert-defined-indirect = $(if $($1),,$(error Variable '$1' must be defined))
$(foreach I,$(IDS), $(call assert-defined-indirect,$I))

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:%.cpp=%.o)

$(OBJS) : CXXFLAGS += $(LIBUTIL_DIR)

$(LIBTRIO): $(OBJS)
	$(AR) -csru $@ $(OBJS)

clean:
	rm -f $(LIBTRIO) $(OBJS)

