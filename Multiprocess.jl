# Funktio prosessien lisäämiseksi
function func_addprocs(nPROCS::Int64)
      if nworkers() == 1 && nPROCS > 1
            addprocs(nPROCS);
      elseif nworkers() < nPROCS
            addprocs(nPROCS-nworkers());
      elseif nworkers() > nPROCS
            rmprocs(workers()[nPROCS+1:end]);
      end
      return nothing
end
