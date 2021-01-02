classdef StateSpaceRepresentation
    properties
        A;
        B;
        C;
        D;        
        system;
    end
    
    methods
        function object = StateSpaceRepresentation(A, B, C, D)
            object.A = A;
            object.B = B;
            object.C = C;
            object.D = D;            
            object.system = ss(object.A, object.B, object.C, object.D);
        end
        
        function setSystemNames(object, stateName, inputName, outputName)
            object.system.StateName = stateName;
            object.system.InputName = inputName;
            object.system.OutputName = outputName;
        end
        
        function response = getSystem(object)
            response = ss(object.A, object.B, object.C, object.D);
        end
        
        %Setters
        function object = set.A(object,A)
            object.A = A; 
            object.system = ss(object.A, object.B, object.C, object.D);
        end
        function object = set.B(object,B)
            object.B = B; 
            object.system = ss(object.A, object.B, object.C, object.D);
        end
        function object = set.C(object,C)
            object.C = C;
            object.system = ss(object.A, object.B, object.C, object.D);
        end
        function object = set.D(object,D)
            object.D = D;
            object.system = ss(object.A, object.B, object.C, object.D);
        end
        %Getters
        function response = get.A(object)
            response = object.A; 
        end
        function response = get.B(object)
            response = object.B; 
        end
        function response = get.C(object)
            response = object.C; 
        end
        function response = get.D(object)
            response = object.D; 
        end
        function response = get.system(object)
            response = object.system; 
        end
    end
end

