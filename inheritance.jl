module inheritance

export Inheritable

abstract Inheritable
 
+(a::Inheritable, b::Inheritable) = (function (a::Inheritable, b::Inheritable)
    properties = Dict{String, Any}()
 
    for property in names(a)
        propertyName = string(property)
 
        if (!haskey(properties, propertyName))
            properties[propertyName] = (propertyName, string(fieldtype(a, property)))
        else
            (fieldName, fieldType) = properties[propertyName]
            properties[propertyName] = (fieldName, "Any")
        end
    end
 
    for property in names(b)
        propertyName = string(property)
 
        if (!haskey(properties, propertyName))
            properties[propertyName] = (propertyName, string(fieldtype(b, property)))
        else
            (fieldName, fieldType) = properties[propertyName]
            properties[propertyName] = (fieldName, "Any")
        end
    end
 
    fieldCode = ""
 
    for property in values(properties)
        (fieldName, fieldType) = property
        fieldCode = fieldCode * fieldName * "::" * fieldType * "\n"
    end
 
    randomTypeName = "An" * randstring(16);
 
    typeCode = "type " * randomTypeName * " <: Inheritable " * fieldCode * " function " * randomTypeName * "() return new () end end"
 
    eval(parse(typeCode))
 
    randomTypeName = symbol(randomTypeName)
 
    c = @eval begin
        $randomTypeName()
    end
 
    for property in names(a)
        try
            c.(property) = a.(property)
        catch
            c
        end
    end
 
    for property in names(b)
        try
            c.(property) = b.(property)
        catch
            c
        end
    end
 
    return c
end)(a, b)

end