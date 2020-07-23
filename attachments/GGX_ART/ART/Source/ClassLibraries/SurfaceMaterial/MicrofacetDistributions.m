/* ===========================================================================

 Copyright (c) 1996-2019 The ART Development Team
 -------------------------------------------

 For a comprehensive list of the members of the development team, and a
 description of their respective contributions, see the file
 "ART_DeveloperList.txt" that is distributed with the libraries.

 This file is part of the Advanced Rendering Toolkit (ART) libraries.

 ART is free software: you can redistribute it and/or modify it under the
 terms of the GNU General Public License as published by the Free Software
 Foundation, either version 3 of the License, or (at your option) any
 later version.

 ART is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 for more details.

 You should have received a copy of the GNU General Public License
 along with ART.  If not, see <http://www.gnu.org/licenses/>.

 =========================================================================== */

#define ART_MODULE_NAME     MicrofacetDistributions

#import "MicrofacetDistributions.h"
#import "SurfaceMaterialMacros.h"

ART_MODULE_INITIALISATION_FUNCTION
(
    [ ArnBlinnMicrofacetDistribution registerWithRuntime ];
    [ ArnGGXMicrofacetDistribution registerWithRuntime ];
)

ART_NO_MODULE_SHUTDOWN_FUNCTION_NECESSARY

// similar to macros in MicrofacetMacros, just with dereferencing
#define MICROFACET_NORMAL_WORLDSPACE ARDIRECTIONCOSINE_VECTOR(*microfacetNormalDirection)
#define MICROFACET_COSINE_WORLDSPACE ARDIRECTIONCOSINE_COSINE(*microfacetNormalDirection)

#define EXPRESSION_FOR_BETA  ((ArNode <ArpDoubleValues> *) ARNUNARY_SUBNODE)

/* ===========================================================================
    'ArnBlinnMicrofacetDistribution'
=========================================================================== */

@implementation ArnBlinnMicrofacetDistribution

ARPCONCRETECLASS_DEFAULT_IMPLEMENTATION(ArnBlinnMicrofacetDistribution)

- init
        : (ArNode *) expressionForBeta
{
    ART_ERRORHANDLING_MANDATORY_ARPROTOCOL_CHECK(
        expressionForBeta,
        ArpDoubleValues
        );
    
    self =
        [ super init
            :   HARD_NODE_REFERENCE(expressionForBeta)
            ];
    
    return self;
}

double  blinn_exponentFromBeta(
        const double  beta
        )
{
    return -log( 2.0 ) / log( cos( M_MAX( 0.0001, beta )));
}

- (double) _sampleProbability
        : (const ArDirectionCosine *) microfacetNormalDirection
        : (const double) exponent
{
    return ( exponent + 1.0 )
        * pow( ARDIRECTIONCOSINE_COSINE(*microfacetNormalDirection), exponent )
        * MATH_1_DIV_2_PI;
}

- (double) _distribution
        : (const ArDirectionCosine *) microfacetNormalDirection
        : (const double) exponent
{
    // the distribution function uses exp+2 even though generation uses exp+1, this is intentional
    return ( exponent + 2.0 )
        * pow( ARDIRECTIONCOSINE_COSINE(*microfacetNormalDirection), exponent )
        * MATH_1_DIV_2_PI;
}

- (BOOL) sampleDistribution
        : (      ArcIntersection *)               incomingDirectionAndLocation
        : (      ArMicrofacetGenerationContext *) context
        : (      ArDirectionCosine *)             sampledNormalDirection
        : (      double *)                        sampleProbability
        : (      double *)                        distributionValue
{
    double  u1, u2;
    
    [ RANDOM_GENERATOR getValuesFromNewSequences: &u1 :&u2 ];

    double  beta;

    [ EXPRESSION_FOR_BETA getDoubleValue
        :   incomingDirectionAndLocation
        : & beta
        ];

    double  exponent = blinn_exponentFromBeta( beta );
    
    double  cosThetaH = pow( u1, 1.0 / ( exponent + 1.0 ));
    ASSERT_UNIT_RANGE_DOUBLE(cosThetaH);

    double sinThetaH = sqrt( 1.0 - M_SQR(cosThetaH) );
    double phi = u2 * MATH_2_MUL_PI;

    Vec3D localH = VEC3D(
        sinThetaH * cos( phi ),
        sinThetaH * sin( phi ),
        cosThetaH
    );
    
    // transform the microfacet normal from local space to world space using
    // the normal of the macrosurface
    Trafo3D     local2world;
    trafo3d_v_local2world_from_worldspace_normal_t(
        & SURFACE_NORMAL_WORLDSPACE,
        & local2world
        );
    vec3d_v_trafo3d_v(
        & localH,
        & local2world,
        & ARDIRECTIONCOSINE_VECTOR(*sampledNormalDirection)
        );
    ARDIRECTIONCOSINE_COSINE(*sampledNormalDirection) = ZC(localH);
    
    // compute the pdf for the microfacet normal from Blinn distribution
    *sampleProbability =
        [ self _sampleProbability
            :  sampledNormalDirection
            :  exponent
            ];
    
    // compute the actual value of the distribution function
    *distributionValue =
        [ self _distribution
            :   sampledNormalDirection
            :   exponent
            ];
    
    return YES;
}

- (BOOL) evaluateDistribution
        : (      ArcIntersection *)               incomingDirectionAndLocation
        : (const ArDirectionCosine *)             microfacetNormalDirection
        : (      ArMicrofacetGenerationContext *) context
        : (      double *)                        sampleProbability
        : (      double *)                        distributionValue
{
    if(ARDIRECTIONCOSINE_COSINE(*microfacetNormalDirection) < 0.0)
        return NO;

    double  beta;

    [ EXPRESSION_FOR_BETA getDoubleValue
        :   incomingDirectionAndLocation
        : & beta
        ];

    double  exponent = blinn_exponentFromBeta( beta );
    
    if(sampleProbability)
    {
        // compute the pdf for the microfacet normal from Blinn distribution
        *sampleProbability =
            [ self _sampleProbability
                :   microfacetNormalDirection
                :   exponent
                ];
    }
    
    // compute the actual value of the distribution function
    *distributionValue =
        [ self _distribution
            :   microfacetNormalDirection
            :   exponent
            ];
    
    return YES;
}

- (double) geometricAttenuation
        : (      ArcIntersection *)               incomingDirectionAndLocation
        : (const ArDirectionCosine *)             microfacetNormalDirection
        : (const ArDirectionCosine *)             outgoingDirection
{
    double incomingDirectionMicrofacetCosine =
        -vec3d_vv_dot(
            & INCOMING_VECTOR_WORLDSPACE,
            & MICROFACET_NORMAL_WORLDSPACE
            );
    double outgoingDirectionMicrofacetCosine =
        vec3d_vv_dot(
            & OUTGOING_VECTOR_WORLDSPACE,
            & MICROFACET_NORMAL_WORLDSPACE
            );
    
    // if either of the directions doesn't end up on the same side of the microfacet
    // as they are of the macrosurface, the geometric attenuation is defined as zero
    if(  incomingDirectionMicrofacetCosine * INCOMING_COSINE_WORLDSPACE <= 0.0
      || outgoingDirectionMicrofacetCosine * OUTGOING_COSINE_WORLDSPACE <= 0.0
      )
        return 0.0;

    return M_MIN3(   1.0,
                   ( 2.0 * M_ABS(MICROFACET_COSINE_WORLDSPACE * OUTGOING_COSINE_WORLDSPACE / outgoingDirectionMicrofacetCosine) ),
                   ( 2.0 * M_ABS(MICROFACET_COSINE_WORLDSPACE * INCOMING_COSINE_WORLDSPACE / incomingDirectionMicrofacetCosine) )
                 );
}

@end



#define EXPRESSION_FOR_ALPHA  ((ArNode <ArpDoubleValues> *) ARNUNARY_SUBNODE)

/* ===========================================================================
    'ArnGGXMicrofacetDistribution'
=========================================================================== */

@implementation ArnGGXMicrofacetDistribution

ARPCONCRETECLASS_DEFAULT_IMPLEMENTATION(ArnGGXMicrofacetDistribution)

- init
        : (ArNode *) expressionForAlpha
{
    ART_ERRORHANDLING_MANDATORY_ARPROTOCOL_CHECK(
        expressionForAlpha,
        ArpDoubleValues
        );
    
    self =
        [ super init
            :   HARD_NODE_REFERENCE(expressionForAlpha)
            ];
    
    return self;
}

- (double) _smith_G1
        : (const Vec3D)                 direction
        : (double)                      cosineDirection
        : (const ArDirectionCosine *)   microfacetNormalDirection
        : (double)                      alpha
{
    // get tan_theta_v first - angle between direction v and macrosurface normal n
    double cos_theta_v = cosineDirection;
    double tan_theta_v = sqrt(1.0 - M_SQR(cos_theta_v))/cos_theta_v;

    // get tan of theta_v ^2 * alpha^2
    double tan_theta_alpha_2 = M_SQR(alpha)*M_SQR(tan_theta_v);

    // the rest of the equation
    return 2.0 / (1.0 + sqrt(1.0 + tan_theta_alpha_2));
}

- (double) _distribution
        : (const ArDirectionCosine *) microfacetNormalDirection
        : (const double) alpha
{
    // get tan_theta_m first - angle between microsurface normal m and macrosurface normal n
    double cos_theta_m = ARDIRECTIONCOSINE_COSINE(*microfacetNormalDirection);
    double tan_theta_m = sqrt(1.0 - M_SQR(cos_theta_m))/cos_theta_m;

    double alpha_2 = M_SQR(alpha);
    // rest of the equation
    double result = alpha_2 / 
        (M_PI * pow(cos_theta_m,4.0) * M_SQR(alpha_2 + M_SQR(tan_theta_m)));

    // Prevent potential numerical issues in other stages of the model
    return (result * cos_theta_m > 1e-20f) ? result : 0.0;
}

- (double) _sampleProbability
        : (const ArDirectionCosine *) microfacetNormalDirection
        : (const double) alpha
{
    double result = 
        [ self _distribution
            :   microfacetNormalDirection
            :   alpha
            ];

    return result * ARDIRECTIONCOSINE_COSINE(*microfacetNormalDirection);
}

- (BOOL) sampleDistribution
        : (      ArcIntersection *)               incomingDirectionAndLocation
        : (      ArMicrofacetGenerationContext *) context
        : (      ArDirectionCosine *)             sampledNormalDirection
        : (      double *)                        sampleProbability
        : (      double *)                        distributionValue
{
    double  u1, u2;
    
    [ RANDOM_GENERATOR getValuesFromNewSequences: &u1 :&u2 ];

    double  alpha;

    [ EXPRESSION_FOR_ALPHA getDoubleValue
        :   incomingDirectionAndLocation
        : & alpha
        ];

    
    // azimuth
    double phi = (2.0 * M_PI) * u2;

    // polar
    double alpha_2 = M_SQR(alpha);

    double tan_theta_m_2 = alpha_2 * u1 / (1.0 - u1);
    double cos_theta = 1.0/(sqrt(1.0 + tan_theta_m_2));
    double cos_theta_2 = M_SQR(cos_theta);

    // Compute probability density of the sampled position
    double temp = 1.0 + tan_theta_m_2 / alpha_2;
    double cos_theta_3 = M_MAX(cos_theta_2 * cos_theta, 1e-20f);
    double sin_theta = sqrt(1.0 - cos_theta_2);

    Vec3D localH = VEC3D(
        sin_theta * cos( phi ),
        sin_theta * sin( phi ),
        cos_theta
    );
    
    // transform the microfacet normal from local space to world space using
    // the normal of the macrosurface
    Trafo3D     local2world;
    trafo3d_v_local2world_from_worldspace_normal_t(
        & SURFACE_NORMAL_WORLDSPACE,
        & local2world
        );
    vec3d_v_trafo3d_v(
        & localH,
        & local2world,
        & ARDIRECTIONCOSINE_VECTOR(*sampledNormalDirection)
        );
    ARDIRECTIONCOSINE_COSINE(*sampledNormalDirection) = ZC(localH);
    
    // compute the pdf for the microfacet normal from Blinn distribution
    *sampleProbability = 1.0/ (M_PI * alpha_2 * cos_theta_3 * M_SQR(temp));

    
    // compute the actual value of the distribution function
    *distributionValue =
        [ self _distribution
            :   sampledNormalDirection
            :   alpha
            ];
    
    return YES;
}

- (BOOL) evaluateDistribution
        : (      ArcIntersection *)               incomingDirectionAndLocation
        : (const ArDirectionCosine *)             microfacetNormalDirection
        : (      ArMicrofacetGenerationContext *) context
        : (      double *)                        sampleProbability
        : (      double *)                        distributionValue
{
    if(ARDIRECTIONCOSINE_COSINE(*microfacetNormalDirection) < 0.0)
        return NO;

    double  alpha;

    [ EXPRESSION_FOR_ALPHA getDoubleValue
        :   incomingDirectionAndLocation
        : & alpha
        ];
    
    if(sampleProbability)
    {
        // compute the pdf for the microfacet normal from Blinn distribution
        *sampleProbability =
            [ self _sampleProbability
                :   microfacetNormalDirection
                :   alpha
                ];
    }
    
    // compute the actual value of the distribution function
    *distributionValue =
        [ self _distribution
            :   microfacetNormalDirection
            :   alpha
            ];
    
    return YES;
}

- (double) geometricAttenuation
        : (      ArcIntersection *)               incomingDirectionAndLocation
        : (const ArDirectionCosine *)             microfacetNormalDirection
        : (const ArDirectionCosine *)             outgoingDirection
{
    double incomingDirectionMicrofacetCosine =
        -vec3d_vv_dot(
            & INCOMING_VECTOR_WORLDSPACE,
            & MICROFACET_NORMAL_WORLDSPACE
            );
    double outgoingDirectionMicrofacetCosine =
        vec3d_vv_dot(
            & OUTGOING_VECTOR_WORLDSPACE,
            & MICROFACET_NORMAL_WORLDSPACE
            );
    
    // if either of the directions doesn't end up on the same side of the microfacet
    // as they are of the macrosurface, the geometric attenuation is defined as zero
    if(  incomingDirectionMicrofacetCosine * INCOMING_COSINE_WORLDSPACE <= 0.0
      || outgoingDirectionMicrofacetCosine * OUTGOING_COSINE_WORLDSPACE <= 0.0
      )
        return 0.0;

    double  alpha;

    [ EXPRESSION_FOR_ALPHA getDoubleValue
        :   incomingDirectionAndLocation
        : & alpha
        ];

    double incomingSmith =
            [ self _smith_G1
                :   INCOMING_VECTOR_WORLDSPACE
                :   INCOMING_COSINE_WORLDSPACE
                :   microfacetNormalDirection
                :   alpha
                ];

    double outgoingSmith =
            [ self _smith_G1
                :   OUTGOING_VECTOR_WORLDSPACE
                :   OUTGOING_COSINE_WORLDSPACE
                :   microfacetNormalDirection
                :   alpha
                ];

    return incomingSmith * outgoingSmith;
}

@end
