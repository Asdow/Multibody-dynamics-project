# Import packages that are needed
import StaticArrays;
import GLVisualize, GLAbstraction, GeometryTypes, Colors;
@everywhere begin
      # Importataan moduulit worker-prosesseille
      if nprocs() > 1
            import StaticArrays;
            import GLVisualize, GLAbstraction, GeometryTypes, Colors;
      end
      global const sa = StaticArrays;
      global const gl = GLVisualize;
      global const gla = GLAbstraction;
      global const gt = GeometryTypes;
      global const col = Colors;
end
