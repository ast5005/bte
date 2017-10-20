
module inheritance_test_types
using inheritance
export A, B
type A <: Inheritable
    whoAmI::Function
    uniqueFunctionA::Function
 
    function A()
        instance = new()
 
        instance.whoAmI = function ()
            println("I am object A")
        end
 
        instance.uniqueFunctionA = function ()
            println("Function unique to A")
        end
 
        return instance
    end
end
 
type B <: Inheritable
    whoAmI::Function
    uniqueFunctionB::Function
 
    function B()
        instance = new()
 
        instance.whoAmI = function ()
            println("I am object B")
        end
 
        instance.uniqueFunctionB = function ()
            println("Function unique to B")
        end
 
        return A() + instance
    end
end

end