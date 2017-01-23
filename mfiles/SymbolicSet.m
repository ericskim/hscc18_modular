classdef SymbolicSet < handle 
% two matlab commands to load grid points from a file stored by
% write_to_file(const SymbolicSet& set, const BDD& bdd, const std::string& filename, char mode='B')
% 
% Please have a look in ./manual/manual.pdf for usage details.
%
  properties (SetAccess=private)
    filename  % name of file which contains the SymbolicSet
    h_ss      % handle to SymbolicSet
    h_bdd     % handle to BDD
    h_cudd    % handle to Cudd
  end
  methods 
    function obj=SymbolicSet(filename)
    % the constructor opens the file 
      if(isstr(filename))
        obj.filename=filename;
      else
        error('filname is not a string');
      end

      [h_ss h_bdd h_cudd]=mexSymbolicSet('init',filename);
      obj.filename=filename;
      obj.h_ss = h_ss;
      obj.h_bdd = h_bdd;
      obj.h_cudd = h_cudd;
      if(isempty(h_ss))
        error(['SymbolicSet: could not read SymbolicSet from ', filename])
      end
    end
    function disp(obj)
      disp(['Matlab object to access the points cointained in the SymbolicSet stored in ', obj.filename])
      disp(' ')
    end
    function delete(obj)
      if(~isempty(obj.h_ss))
        mexSymbolicSet('delete', obj.h_ss,obj.h_bdd, obj.h_cudd);
      end
    end


    %function points=get.points(obj)
    %% read the grid points from file
    %% optional projection indices
    %% get.point('project',[1 2 3]) loads the grid points in the symbolic set
    %% projected onto the dimensions [1 2 3]

    %  if(obj.project)
    %    project=obj.project-1;
    %    points=mexSymbolicSet(obj.filename,'gridpoints',project);
    %    obj.points=points;
    %    if(points==0)
    %      error(['there are not points in the Symbolic Set ',obj.filename]);
    %    end
    %  else
    %    points=mexSymbolicSet(obj.filename,'gridpoints');
    %    obj.points=points;
    %    if(points==0)
    %      error(['there are not points in the Symbolic Set ',obj.filename]);
    %    end
    %  end
    %end
  end
end
