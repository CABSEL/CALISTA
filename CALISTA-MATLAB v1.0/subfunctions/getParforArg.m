function arg = getParforArg()
p = gcp('nocreate'); % get the pool object if it exists, but never open a pool
if isempty(p)
  arg = 0;
else
  arg = p.NumWorkers;
end
end
