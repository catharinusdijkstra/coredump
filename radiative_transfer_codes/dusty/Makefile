# Compiler
FC = gfortran

# Compiler flags
FFLAGS = -std=legacy

# Source files
SOURCES = dusty.f input_output.f radiative_transfer.f radiatively_driven_winds_dynamics.f \
          slab_geometry.f radiative_transfer_results_analysis.f \
          math_and_number_to_string_conversion.f

# Executable name
TARGET = dusty.exe

# Build rule
$(TARGET): $(SOURCES)
	$(FC) $(FFLAGS) -o $@ $^

# Clean rule
clean:
	rm -f $(TARGET) *.o
