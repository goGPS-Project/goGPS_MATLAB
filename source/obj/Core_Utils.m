%   CLASS Core_Utils
% =========================================================================
%
% DESCRIPTION
%   Class to manages utilities
%
% EXAMPLE
%   % set of static utilities
%   Core_Utils.diffAndPred
%
% FOR A LIST OF CONSTANTs and METHODS use doc Core_UI

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
%  Written by:        Andrea Gatti, Giulio Tagliaferro
%  Contributors:      Andrea Gatti, Giulio Tagliaferro, ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

classdef Core_Utils < handle
    properties (Constant)
        V_LIGHT = 299792458;                % Velocity of light in the void [m/s]
    end
    
    properties (Constant, Access = private)
        ISO_DEF = ['Country','Alpha-2 code','Alpha-3 code','Numeric code','Latitude (average)','Longitude (average)'];
        ISO3166 = { ...
            'AF', 'AFG', 'Afghanistan', '4', '33';
            'AL', 'ALB', 'Albania', '8', '41';
            'DZ', 'DZA', 'Algeria', '12', '28';
            'AS', 'ASM', 'American Samoa', '16', '-14.3333';
            'AD', 'AND', 'Andorra', '20', '42.5';
            'AO', 'AGO', 'Angola', '24', '-12.5';
            'AI', 'AIA', 'Anguilla', '660', '18.25';
            'AQ', 'ATA', 'Antarctica', '10', '-90';
            'AG', 'ATG', 'Antigua and Barbuda', '28', '17.05';
            'AR', 'ARG', 'Argentina', '32', '-34';
            'AM', 'ARM', 'Armenia', '51', '40';
            'AW', 'ABW', 'Aruba', '533', '12.5';
            'AU', 'AUS', 'Australia', '36', '-27';
            'AT', 'AUT', 'Austria', '40', '47.3333';
            'AZ', 'AZE', 'Azerbaijan', '31', '40.5';
            'BS', 'BHS', 'Bahamas', '44', '24.25';
            'BH', 'BHR', 'Bahrain', '48', '26';
            'BD', 'BGD', 'Bangladesh', '50', '24';
            'BB', 'BRB', 'Barbados', '52', '13.1667';
            'BY', 'BLR', 'Belarus', '112', '53';
            'BE', 'BEL', 'Belgium', '56', '50.8333';
            'BZ', 'BLZ', 'Belize', '84', '17.25';
            'BJ', 'BEN', 'Benin', '204', '9.5';
            'BM', 'BMU', 'Bermuda', '60', '32.3333';
            'BT', 'BTN', 'Bhutan', '64', '27.5';
            'BO', 'BOL', 'Bolivia, Plurinational State of', '68', '-17';
            'BO', 'BOL', 'Bolivia', '68', '-17';
            'BA', 'BIH', 'Bosnia and Herzegovina', '70', '44';
            'BW', 'BWA', 'Botswana', '72', '-22';
            'BV', 'BVT', 'Bouvet Island', '74', '-54.4333';
            'BR', 'BRA', 'Brazil', '76', '-10';
            'IO', 'IOT', 'British Indian Ocean Territory', '86', '-6';
            'BN', 'BRN', 'Brunei Darussalam', '96', '4.5';
            'BN', 'BRN', 'Brunei', '96', '4.5';
            'BG', 'BGR', 'Bulgaria', '100', '43';
            'BF', 'BFA', 'Burkina Faso', '854', '13';
            'BI', 'BDI', 'Burundi', '108', '-3.5';
            'KH', 'KHM', 'Cambodia', '116', '13';
            'CM', 'CMR', 'Cameroon', '120', '6';
            'CA', 'CAN', 'Canada', '124', '60';
            'CV', 'CPV', 'Cape Verde', '132', '16';
            'KY', 'CYM', 'Cayman Islands', '136', '19.5';
            'CF', 'CAF', 'Central African Republic', '140', '7';
            'TD', 'TCD', 'Chad', '148', '15';
            'CL', 'CHL', 'Chile', '152', '-30';
            'CN', 'CHN', 'China', '156', '35';
            'CX', 'CXR', 'Christmas Island', '162', '-10.5';
            'CC', 'CCK', 'Cocos (Keeling) Islands', '166', '-12.5';
            'CO', 'COL', 'Colombia', '170', '4';
            'KM', 'COM', 'Comoros', '174', '-12.1667';
            'CG', 'COG', 'Congo', '178', '-1';
            'CD', 'COD', 'Congo, the Democratic Republic of the', '180', '0';
            'CK', 'COK', 'Cook Islands', '184', '-21.2333';
            'CR', 'CRI', 'Costa Rica', '188', '10';
            'CI', 'CIV', 'Côte d''Ivoire', '384', '8';
            'CI', 'CIV', 'Ivory Coast', '384', '8';
            'HR', 'HRV', 'Croatia', '191', '45.1667';
            'CU', 'CUB', 'Cuba', '192', '21.5';
            'CY', 'CYP', 'Cyprus', '196', '35';
            'CZ', 'CZE', 'Czech Republic', '203', '49.75';
            'DK', 'DNK', 'Denmark', '208', '56';
            'DJ', 'DJI', 'Djibouti', '262', '11.5';
            'DM', 'DMA', 'Dominica', '212', '15.4167';
            'DO', 'DOM', 'Dominican Republic', '214', '19';
            'EC', 'ECU', 'Ecuador', '218', '-2';
            'EG', 'EGY', 'Egypt', '818', '27';
            'SV', 'SLV', 'El Salvador', '222', '13.8333';
            'GQ', 'GNQ', 'Equatorial Guinea', '226', '2';
            'ER', 'ERI', 'Eritrea', '232', '15';
            'EE', 'EST', 'Estonia', '233', '59';
            'ET', 'ETH', 'Ethiopia', '231', '8';
            'FK', 'FLK', 'Falkland Islands (Malvinas)', '238', '-51.75';
            'FO', 'FRO', 'Faroe Islands', '234', '62';
            'FJ', 'FJI', 'Fiji', '242', '-18';
            'FI', 'FIN', 'Finland', '246', '64';
            'FR', 'FRA', 'France', '250', '46';
            'GF', 'GUF', 'French Guiana', '254', '4';
            'PF', 'PYF', 'French Polynesia', '258', '-15';
            'TF', 'ATF', 'French Southern Territories', '260', '-43';
            'GA', 'GAB', 'Gabon', '266', '-1';
            'GM', 'GMB', 'Gambia', '270', '13.4667';
            'GE', 'GEO', 'Georgia', '268', '42';
            'DE', 'DEU', 'Germany', '276', '51';
            'GH', 'GHA', 'Ghana', '288', '8';
            'GI', 'GIB', 'Gibraltar', '292', '36.1833';
            'GR', 'GRC', 'Greece', '300', '39';
            'GL', 'GRL', 'Greenland', '304', '72';
            'GD', 'GRD', 'Grenada', '308', '12.1167';
            'GP', 'GLP', 'Guadeloupe', '312', '16.25';
            'GU', 'GUM', 'Guam', '316', '13.4667';
            'GT', 'GTM', 'Guatemala', '320', '15.5';
            'GG', 'GGY', 'Guernsey', '831', '49.5';
            'GN', 'GIN', 'Guinea', '324', '11';
            'GW', 'GNB', 'Guinea-Bissau', '624', '12';
            'GY', 'GUY', 'Guyana', '328', '5';
            'HT', 'HTI', 'Haiti', '332', '19';
            'HM', 'HMD', 'Heard Island and McDonald Islands', '334', '-53.1';
            'VA', 'VAT', 'Holy See (Vatican City State)', '336', '41.9';
            'HN', 'HND', 'Honduras', '340', '15';
            'HK', 'HKG', 'Hong Kong', '344', '22.25';
            'HU', 'HUN', 'Hungary', '348', '47';
            'IS', 'ISL', 'Iceland', '352', '65';
            'IN', 'IND', 'India', '356', '20';
            'ID', 'IDN', 'Indonesia', '360', '-5';
            'IR', 'IRN', 'Iran, Islamic Republic of', '364', '32';
            'IQ', 'IRQ', 'Iraq', '368', '33';
            'IE', 'IRL', 'Ireland', '372', '53';
            'IM', 'IMN', 'Isle of Man', '833', '54.23';
            'IL', 'ISR', 'Israel', '376', '31.5';
            'IT', 'ITA', 'Italy', '380', '42.8333';
            'JM', 'JAM', 'Jamaica', '388', '18.25';
            'JP', 'JPN', 'Japan', '392', '36';
            'JE', 'JEY', 'Jersey', '832', '49.21';
            'JO', 'JOR', 'Jordan', '400', '31';
            'KZ', 'KAZ', 'Kazakhstan', '398', '48';
            'KE', 'KEN', 'Kenya', '404', '1';
            'KI', 'KIR', 'Kiribati', '296', '1.4167';
            'KP', 'PRK', 'Korea, Democratic People''s Republic of', '408', '40';
            'KR', 'KOR', 'Korea, Republic of', '410', '37';
            'KR', 'KOR', 'South Korea', '410', '37';
            'KW', 'KWT', 'Kuwait', '414', '29.3375';
            'KG', 'KGZ', 'Kyrgyzstan', '417', '41';
            'LA', 'LAO', 'Lao People''s Democratic Republic', '418', '18';
            'LV', 'LVA', 'Latvia', '428', '57';
            'LB', 'LBN', 'Lebanon', '422', '33.8333';
            'LS', 'LSO', 'Lesotho', '426', '-29.5';
            'LR', 'LBR', 'Liberia', '430', '6.5';
            'LY', 'LBY', 'Libyan Arab Jamahiriya', '434', '25';
            'LY', 'LBY', 'Libya', '434', '25';
            'LI', 'LIE', 'Liechtenstein', '438', '47.1667';
            'LT', 'LTU', 'Lithuania', '440', '56';
            'LU', 'LUX', 'Luxembourg', '442', '49.75';
            'MO', 'MAC', 'Macao', '446', '22.1667';
            'MK', 'MKD', 'Macedonia, the former Yugoslav Republic of', '807', '41.8333';
            'MG', 'MDG', 'Madagascar', '450', '-20';
            'MW', 'MWI', 'Malawi', '454', '-13.5';
            'MY', 'MYS', 'Malaysia', '458', '2.5';
            'MV', 'MDV', 'Maldives', '462', '3.25';
            'ML', 'MLI', 'Mali', '466', '17';
            'MT', 'MLT', 'Malta', '470', '35.8333';
            'MH', 'MHL', 'Marshall Islands', '584', '9';
            'MQ', 'MTQ', 'Martinique', '474', '14.6667';
            'MR', 'MRT', 'Mauritania', '478', '20';
            'MU', 'MUS', 'Mauritius', '480', '-20.2833';
            'YT', 'MYT', 'Mayotte', '175', '-12.8333';
            'MX', 'MEX', 'Mexico', '484', '23';
            'FM', 'FSM', 'Micronesia, Federated States of', '583', '6.9167';
            'MD', 'MDA', 'Moldova, Republic of', '498', '47';
            'MC', 'MCO', 'Monaco', '492', '43.7333';
            'MN', 'MNG', 'Mongolia', '496', '46';
            'ME', 'MNE', 'Montenegro', '499', '42';
            'MS', 'MSR', 'Montserrat', '500', '16.75';
            'MA', 'MAR', 'Morocco', '504', '32';
            'MZ', 'MOZ', 'Mozambique', '508', '-18.25';
            'MM', 'MMR', 'Myanmar', '104', '22';
            'MM', 'MMR', 'Burma', '104', '22';
            'NA', 'NAM', 'Namibia', '516', '-22';
            'NR', 'NRU', 'Nauru', '520', '-0.5333';
            'NP', 'NPL', 'Nepal', '524', '28';
            'NL', 'NLD', 'Netherlands', '528', '52.5';
            'AN', 'ANT', 'Netherlands Antilles', '530', '12.25';
            'NC', 'NCL', 'New Caledonia', '540', '-21.5';
            'NZ', 'NZL', 'New Zealand', '554', '-41';
            'NI', 'NIC', 'Nicaragua', '558', '13';
            'NE', 'NER', 'Niger', '562', '16';
            'NG', 'NGA', 'Nigeria', '566', '10';
            'NU', 'NIU', 'Niue', '570', '-19.0333';
            'NF', 'NFK', 'Norfolk Island', '574', '-29.0333';
            'MP', 'MNP', 'Northern Mariana Islands', '580', '15.2';
            'NO', 'NOR', 'Norway', '578', '62';
            'OM', 'OMN', 'Oman', '512', '21';
            'PK', 'PAK', 'Pakistan', '586', '30';
            'PW', 'PLW', 'Palau', '585', '7.5';
            'PS', 'PSE', 'Palestinian Territory, Occupied', '275', '32';
            'PA', 'PAN', 'Panama', '591', '9';
            'PG', 'PNG', 'Papua New Guinea', '598', '-6';
            'PY', 'PRY', 'Paraguay', '600', '-23';
            'PE', 'PER', 'Peru', '604', '-10';
            'PH', 'PHL', 'Philippines', '608', '13';
            'PN', 'PCN', 'Pitcairn', '612', '-24.7';
            'PL', 'POL', 'Poland', '616', '52';
            'PT', 'PRT', 'Portugal', '620', '39.5';
            'PR', 'PRI', 'Puerto Rico', '630', '18.25';
            'QA', 'QAT', 'Qatar', '634', '25.5';
            'RE', 'REU', 'Réunion', '638', '-21.1';
            'RO', 'ROU', 'Romania', '642', '46';
            'RU', 'RUS', 'Russian Federation', '643', '60';
            'RU', 'RUS', 'Russia', '643', '60';
            'RW', 'RWA', 'Rwanda', '646', '-2';
            'SH', 'SHN', 'Saint Helena, Ascension and Tristan da Cunha', '654', '-15.9333';
            'KN', 'KNA', 'Saint Kitts and Nevis', '659', '17.3333';
            'LC', 'LCA', 'Saint Lucia', '662', '13.8833';
            'PM', 'SPM', 'Saint Pierre and Miquelon', '666', '46.8333';
            'VC', 'VCT', 'Saint Vincent and the Grenadines', '670', '13.25';
            'VC', 'VCT', 'Saint Vincent & the Grenadines', '670', '13.25';
            'VC', 'VCT', 'St. Vincent and the Grenadines', '670', '13.25';
            'WS', 'WSM', 'Samoa', '882', '-13.5833';
            'SM', 'SMR', 'San Marino', '674', '43.7667';
            'ST', 'STP', 'Sao Tome and Principe', '678', '1';
            'SA', 'SAU', 'Saudi Arabia', '682', '25';
            'SN', 'SEN', 'Senegal', '686', '14';
            'RS', 'SRB', 'Serbia', '688', '44';
            'SC', 'SYC', 'Seychelles', '690', '-4.5833';
            'SL', 'SLE', 'Sierra Leone', '694', '8.5';
            'SG', 'SGP', 'Singapore', '702', '1.3667';
            'SK', 'SVK', 'Slovakia', '703', '48.6667';
            'SI', 'SVN', 'Slovenia', '705', '46';
            'SB', 'SLB', 'Solomon Islands', '90', '-8';
            'SO', 'SOM', 'Somalia', '706', '10';
            'ZA', 'ZAF', 'South Africa', '710', '-29';
            'GS', 'SGS', 'South Georgia and the South Sandwich Islands', '239', '-54.5';
            'SS', 'SSD', 'South Sudan', '728', '8';
            'ES', 'ESP', 'Spain', '724', '40';
            'LK', 'LKA', 'Sri Lanka', '144', '7';
            'SD', 'SDN', 'Sudan', '736', '15';
            'SR', 'SUR', 'Suriname', '740', '4';
            'SJ', 'SJM', 'Svalbard and Jan Mayen', '744', '78';
            'SZ', 'SWZ', 'Swaziland', '748', '-26.5';
            'SE', 'SWE', 'Sweden', '752', '62';
            'CH', 'CHE', 'Switzerland', '756', '47';
            'SY', 'SYR', 'Syrian Arab Republic', '760', '35';
            'TW', 'TWN', 'Taiwan, Province of China', '158', '23.5';
            'TW', 'TWN', 'Taiwan', '158', '23.5';
            'TJ', 'TJK', 'Tajikistan', '762', '39';
            'TZ', 'TZA', 'Tanzania, United Republic of', '834', '-6';
            'TH', 'THA', 'Thailand', '764', '15';
            'TL', 'TLS', 'Timor-Leste', '626', '-8.55';
            'TG', 'TGO', 'Togo', '768', '8';
            'TK', 'TKL', 'Tokelau', '772', '-9';
            'TO', 'TON', 'Tonga', '776', '-20';
            'TT', 'TTO', 'Trinidad and Tobago', '780', '11';
            'TN', 'TUN', 'Tunisia', '788', '34';
            'TR', 'TUR', 'Turkey', '792', '39';
            'TM', 'TKM', 'Turkmenistan', '795', '40';
            'TC', 'TCA', 'Turks and Caicos Islands', '796', '21.75';
            'TV', 'TUV', 'Tuvalu', '798', '-8';
            'UG', 'UGA', 'Uganda', '800', '1';
            'UA', 'UKR', 'Ukraine', '804', '49';
            'AE', 'ARE', 'United Arab Emirates', '784', '24';
            'GB', 'GBR', 'United Kingdom', '826', '54';
            'US', 'USA', 'United States', '840', '38';
            'UM', 'UMI', 'United States Minor Outlying Islands', '581', '19.2833';
            'UY', 'URY', 'Uruguay', '858', '-33';
            'UZ', 'UZB', 'Uzbekistan', '860', '41';
            'VU', 'VUT', 'Vanuatu', '548', '-16';
            'VE', 'VEN', 'Venezuela, Bolivarian Republic of', '862', '8';
            'VE', 'VEN', 'Venezuela', '862', '8';
            'VN', 'VNM', 'Viet Nam', '704', '16';
            'VN', 'VNM', 'Vietnam', '704', '16';
            'VG', 'VGB', 'Virgin Islands, British', '92', '18.5';
            'VI', 'VIR', 'Virgin Islands, U.S.', '850', '18.3333';
            'WF', 'WLF', 'Wallis and Futuna', '876', '-13.3';
            'EH', 'ESH', 'Western Sahara', '732', '24.5';
            'YE', 'YEM', 'Yemen', '887', '15';
            'ZM', 'ZMB', 'Zambia', '894', '-15';
            'ZW', 'ZWE', 'Zimbabwe', '716', '-20';
            'OC', 'OCE', 'Ocean', '0', '0';   % fake entry
            'WD', 'WRD', 'World', '-90', '0'; % fake entry
            };
    end
    
    methods (Static)               
        %--------------------------------------------------------------------------
        %% FILTERS AND INTERPOLATORS
        %--------------------------------------------------------------------------
        
        function diff_data = diffAndPred(data, n_order, t_ref, method)
            % compute diff predicting epoch 0 of each arc
            % using interp 1 spline method
            %
            % SYNTAX
            %   Core_Utils.diffAndPred(data, n_order, t_ref, method)
            
            if nargin < 3 || isempty(t_ref)
                t_ref = 1 : size(data,1);
            end
            if nargin < 2 || isempty(n_order)
                n_order = 1;
            end
            if nargin < 4 || isempty(method)
                method = 'pchip';
            end
            if isempty(data)
                diff_data = [];
            else
                diff_data = nan(size(data));
                % Add n_order rows to data
                data = [repmat(data(1,:), n_order, 1); data];
                for s = 1 : size(data, 2)
                    % Get the original good data for column s
                    tmp = data(1 + n_order : end, s);
                    id_ok = ~isnan(tmp);
                    if sum(id_ok) > 2
                        lim = getFlagsLimits(id_ok);
                        % Interpolate data beginning
                        % interpolate the "left" of the first element of an arc
                        % because diff "eat" the first value
                        %if (length(id_ok) > (n_order + 1)) && any(id_ok(1))
                        %    id_est = find(id_ok(lim(1,1):lim(1,2)));
                        %    data(1 : n_order, s) = interp1(t_ref(id_est), tmp(id_est), 1 - n_order : 0, 'spline', 'extrap');
                        %end
                        
                        lim_short = lim(lim(:,2) - lim(:,1) < 2 & lim(:,1) > 1, :);
                        % short arcs cannot be differenciated efficiently
                        for l = 1 : size(lim_short, 1)
                            data(lim_short(l, 1), s) = data(lim_short(l, 1)+1, s);
                        end
                        
                        % differenciate only limits larger than 2
                        lim = lim(lim(:,2) - lim(:,1) > 1, :);
                        for l = 1 : size(lim, 1)
                            id_data = lim(l, 1) : lim(l, 2);
                            id_est = 0 : (n_order - 1);
                            
                            % slower approach with interp1
                            % data(lim(l, 1) + id_est, s) = interp1(t_ref(id_data), tmp(id_data), lim(l, 1) - 1 - fliplr(id_est), 'spline', 'extrap');
                            
                            % faster approach skipping a lot of checks
                            % this is the internal implementation of interp1
                            if strcmp(method, 'zeros')
                                data(lim(l, 1) + id_est, s) = 0;
                            else
                                fun = griddedInterpolant(t_ref(id_data), tmp(id_data), method);
                                data(lim(l, 1) + id_est, s) = fun(lim(l, 1) - 1 - fliplr(id_est));
                            end
                            
                            diff_data(id_data, s) = diff(data(lim(l, 1) : (lim(l, 2) + n_order), s), n_order);
                            % restore data for the next interval
                            data(1 + n_order : end, s) = tmp;
                        end
                    end
                end
                % diff_data = diff(data, n_order); % now it is done arc by arc
            end
        end
        
        
        function [y_out,coeff] = interp1LS(x_in, y_in, degree, x_out, flag_reg)
            % Least squares interpolant of a 1D dataset
            %
            % SYNTAX
            %   y_out = interp1LS(x_in, y_in, degree, x_out)
            
            if nargin < 4 || isempty(x_out)
                x_out = x_in;
            end
            if nargin < 5
                flag_reg = false;
            end
            coeff = nan(degree+1, iif(min(size(y_in)) == 1, 1, size(y_in,2)));
            for c = 1 : iif(min(size(y_in)) == 1, 1, size(y_in,2))
                if size(y_in, 1) == 1
                    y_tmp = y_in';
                else
                    y_tmp = y_in(:, c);
                end
                id_ok = ~isnan(y_tmp(:)) & ~isnan(x_in(:));
                x_tmp = x_in(id_ok);
                y_tmp = y_tmp(id_ok);
                
                n_obs = numel(x_tmp);
                A = zeros(n_obs, degree + 1);
                A(:, 1) = ones(n_obs, 1);
                for d = 1 : degree
                    A(:, d + 1) = x_tmp .^ d;
                end
                
                if (nargin < 4) && numel(x_out) == numel(x_tmp) &&  (sum(x_out(:) - x_tmp(:)) == 0)
                    A2 = A;
                else
                    n_out = numel(x_out);
                    A2 = zeros(n_out, degree + 1);
                    A2(:, 1) = ones(n_out, 1);
                    for d = 1 : degree
                        A2(:, d + 1) = x_out .^ d;
                    end
                end
                
                warning('off')
                if (min(size(y_in)) == 1) && (numel(flag_reg) == 1) && not(flag_reg)
                    %y_out = A2 * ((A' * A) \ (A' * y_tmp(:)));
                    coeff(:,c) = (A \ y_tmp(:));
                    y_out = A2 * coeff(:,c);
                    y_out = reshape(y_out, size(x_out, 1), size(x_out, 2));
                elseif numel(flag_reg) == numel(y_in)
                    if size(flag_reg, 1) == 1
                        v_tmp = flag_reg';
                    else
                        v_tmp = flag_reg(:, c);
                    end
                    % full regularization
                    invQ = spdiags(1 ./ (v_tmp.^2), 0, numel(v_tmp), numel(v_tmp));
                    At_invQ = A' * invQ;
                    coeff(:,c) = full(((At_invQ * A) \ (At_invQ * y_tmp(:))));
                    y_out(:,c) = A2 * coeff(:,c);
                else
                    coeff(:,c) = full(((A' * A + 1e-6 * speye(size(A,2))) \ (A' * y_tmp(:))));
                    y_out(:,c) = A2 * coeff(:,c);
                end
                warning('on')
            end
        end
        
        function [y_out] = interpPlaneLS(x_in, y_in, degree, x_out)
            % Least squares interpolant of a 3D dataset
            % Estimate a plane on x, y + polynomial on z
            % y_in is an array containing m set of data (1 set per column)
            % e.g. degree = 2
            %   f(x,y,z) =  a*x + b*y + c + d*z^2 + e*z
            %
            % SYNTAX
            %   y_out = interpPlane1LS(x_in, y_in, z_degree, x_out)
            
            if nargin < 4
                x_out = x_in;
            end
            
            for c = 1 : iif(min(size(y_in)) == 1, 1, size(y_in,2))
                % try to correct orientation of x
                if size(y_in, 1) == 1 && size(y_tmp, 1) ~= size(x_tmp, 1)
                    y_tmp = y_in';
                else
                    y_tmp = y_in(:, c);
                end
                id_ok = ~isnan(y_tmp(:)) & ~any(isnan(x_in), 2);
                x_tmp = x_in(id_ok, :);
                y_tmp = y_tmp(id_ok);
                
                n_obs = size(x_tmp, 1);
                A = zeros(n_obs, degree + 3);
                A(:, 1) = ones(n_obs, 1);
                A(:, 2) = x_tmp(:, 1); % y
                A(:, 3) = x_tmp(:, 2); % x
                for d = 1 : degree
                    A(:, d + 3) = x_tmp(:, 3) .^ d;
                end
                
                if (nargin < 4) && numel(x_out) == numel(x_tmp) &&  (sum(x_out(:) - x_tmp(:)) == 0)
                    A2 = A;
                else
                    n_out = size(x_out, 1);
                    A2 = zeros(n_out, degree + 3);
                    A2(:, 1) = ones(n_out, 1);
                    A2(:, 2) = x_out(:, 1); % y
                    A2(:, 3) = x_out(:, 2); % x
                    for d = 1 : degree
                        A2(:, d + 3) = x_out(:, 3) .^ d;
                    end
                end
                
                warning('off')
                if min(size(y_in)) == 1
                    y_out = A2 * ((A' * A) \ (A' * y_tmp(:)));
                    y_out = reshape(y_out, size(x_out, 1), size(x_out, 3));
                else
                    y_out(:,c) = A2 * ((A' * A + eye(size(A,2))) \ (A' * y_tmp(:)));
                end
                warning('on')
            end
        end
        
        function data_out = resize2(data, size_out)
            % Simple resize of a 2D dataset
            %
            % SYNTAX 
            %   data_out = resize2(data, size_out)
            if all(size(data) == size_out)
                data_out = data;
            else
                if size(data, 1) == 1
                    data = repmat(data, 2, 1);
                end
                if size(data, 2) == 1
                    data = repmat(data, 1, 2);
                end
                size_in = size(data);
                y_grid = linspace(0 + 0.5/size_in(1), 1 - 0.5/size_in(1), size_in(1));
                x_grid = linspace(0 + 0.5/size_in(2), 1 - 0.5/size_in(2), size_in(2));
                [X, Y] = ndgrid(y_grid, x_grid);
                y_grid_out = linspace(0 + 0.5/size_out(1), 1 - 0.5/size_out(1), size_out(1));
                x_grid_out = linspace(0 + 0.5/size_out(2), 1 - 0.5/size_out(2), size_out(2));
                fInterp = griddedInterpolant(X, Y, data, 'spline');
                [X, Y] = ndgrid(y_grid_out, x_grid_out);
                data_out = fInterp(X, Y);
            end
        end
        
        function [y_out] = interp22nLS(x_in, y_in, degree, x_out)
            % Least squares interpolant of a 3D dataset
            % Estimate a plane on x, y + polynomial on z
            % y_in is an array containing m set of data (1 set per column)
            % e.g. degree = 2
            %   f(x,y,z) =  a*x + b*y + c + d*z^2 + e*z
            %
            % SYNTAX
            %   y_out = interpPlane1LS(x_in, y_in, z_degree, x_out)
            
            if nargin < 4
                x_out = x_in;
            end
            
            for c = 1 : iif(min(size(y_in)) == 1, 1, size(y_in,2))
                % try to correct orientation of x
                if size(y_in, 1) == 1 && size(y_tmp, 1) ~= size(x_tmp, 1)
                    y_tmp = y_in';
                else
                    y_tmp = y_in(:, c);
                end
                id_ok = ~isnan(y_tmp(:)) & ~any(isnan(x_in), 2);
                x_tmp = x_in(id_ok, :);
                y_tmp = y_tmp(id_ok);
                
                n_obs = size(x_tmp, 1);
                A = zeros(n_obs, degree + 6);
                A(:, 1) = ones(n_obs, 1);
                A(:, 2) = x_tmp(:, 1); % y
                A(:, 3) = x_tmp(:, 2); % x
                A(:, 4) = x_tmp(:, 1) .^ 2; % x^2
                A(:, 5) = x_tmp(:, 2) .^ 2; % y^2
                A(:, 6) = x_tmp(:, 1) .* x_tmp(:, 2); % x*y
                for d = 1 : degree
                    A(:, d + 6) = x_tmp(:, 3) .^ d;
                end
                
                if (nargin < 4) && numel(x_out) == numel(x_tmp) &&  (sum(x_out(:) - x_tmp(:)) == 0)
                    A2 = A;
                else
                    n_out = size(x_out, 1);
                    A2 = zeros(n_out, degree + 6);
                    A2(:, 1) = ones(n_out, 1);
                    A2(:, 2) = x_out(:, 1); % y
                    A2(:, 3) = x_out(:, 2); % x
                    A2(:, 4) = x_out(:, 1) .^ 2; % x^2
                    A2(:, 5) = x_out(:, 2) .^ 2; % y^2
                    A2(:, 6) = x_out(:, 1) .* x_out(:, 2); % x*y
                    for d = 1 : degree
                        A2(:, d + 6) = x_out(:, 3) .^ d;
                    end
                end
                
                warning('off')
                if min(size(y_in)) == 1
                    y_out = A2 * ((A' * A) \ (A' * y_tmp(:)));
                    y_out = reshape(y_out, size(x_out, 1), size(x_out, 3));
                else
                    y_out(:,c) = A2 * ((A' * A + eye(size(A,2))) \ (A' * y_tmp(:)));
                end
                warning('on')
            end
        end
        
        function val = linInterpLatLonTime(data, first_lat, dlat, first_lon, dlon, first_t, dt, lat, lon, t)
            % Interpolate values froma data on a gepgraphical grid with multiple epoch
            % data structure:
            %        first dimension : dlat (+) south pole -> north pole
            %        second dimension : dlon (+) west -> east
            %        third dimension : dr (+) time usual direction
            %        NOTE: dlat, dlon,dt do not have to be positive
            %
            % INPUT:
            %      data - the data to be interpolate
            %      fist_lat - value of first lat value (max lat)
            %      dlat - px size lat
            %      first_lon - value of first lon value
            %      dlon - px size lon
            %      first_t - value of first time
            %      dt - px size time
            %      lat - lat at what we want to interpolate
            %      lon - lon at what we ant to interpolate
            %      gps_time - time at what we want to interpolate
            % NOTES 1 - all lat values should have same unit of measure
            %       2 - all lon values should have same unit of measure
            %       3 - all time values should have same unit of measure
            %       4 - the method will interpolate first in the dimesnion with less time
            % IMPORTANT : no double values at the borders should coexist: e.g. -180 180 or 0 360
            [nlat , nlon, nt] = size(data);
            n_in_lat = length(lat);
            n_in_lon = length(lon);
            n_in_t = length(t);
            assert(n_in_lat == n_in_lon);
            [ it, st, ilons, ilone, slon, ilat, slat] = Core_Utils.getIntIdx(data, first_lat, dlat, first_lon, dlon, first_t, dt, lat, lon,t);
            if n_in_lat > n_in_t % time first
                
                it = it*ones(size(ilat));
                % interpolate along time
                % [ 1 2  <= index of the cell at the smae time
                %   3 4]
                idx1 = sub2ind([nlat nlon nt], ilat, ilons, it);
                idx2 = sub2ind([nlat nlon nt], ilat, ilons, it+1);
                vallu = data(idx1).*(1-st) + data(idx2).*st;
                idx1 = sub2ind([nlat nlon nt], ilat   , ilone , it);
                idx2 = sub2ind([nlat nlon nt],ilat   , ilone , it+1);
                valru = data(idx1).*(1-st) + data(idx2).*st;
                idx1 = sub2ind([nlat nlon nt],ilat+1 , ilons , it);
                idx2 = sub2ind([nlat nlon nt],ilat+1 , ilons , it+1);
                valld =  data(idx1).*(1-st) + data(idx2).*st;
                idx1 = sub2ind([nlat nlon nt],ilat+1 , ilone , it);
                idx2 = sub2ind([nlat nlon nt],ilat+1 , ilone , it+1);
                valrd =  data(idx1).*(1-st) + data(idx2).*st;
                
                %interpolate along long
                valu = vallu.*(1-slon) + valru.*slon;
                vald = valld.*(1-slon) + valrd.*slon;
                
                %interpolate along lat
                val = valu.*(1-slat) + vald.*slat;
                
            else %space first % NOTE: consider speed up in case only one time is present, unnecessary operations done                
                % lon first % NOTE: consider speed up in case only one time is present, unnecessary operations done
                % interpolate along lon
                % before up
                % after up
                % before down
                % after down
                if numel(it) == 1
                    it = it*ones(size(ilat));
                end
                idx1 = sub2ind([nlat nlon nt], ilat, ilons(:,1), it);
                idx2 = sub2ind([nlat nlon nt], ilat, ilone(:,1), it);
                valbu = permute(data(idx1).*(1-slon(:,1)) + data(idx2).*slon(:,1),[3 1 2]);
                idx1 = sub2ind([nlat nlon nt], ilat   , ilons(:,1) , min(it+1,size(data,3)));
                idx2 = sub2ind([nlat nlon nt], ilat   , ilone(:,1) , min(it+1,size(data,3)));
                valau = permute(data(idx1).*(1-slon(:,1)) + data(idx2).*slon(:,1),[3 1 2]);
                idx1 = sub2ind([nlat nlon nt], ilat+1 , ilons(:,1) , it);
                idx2 = sub2ind([nlat nlon nt], ilat+1 , ilone(:,1) , it);
                valbd = permute(data(idx1).*(1-slon(:,1)) + data(idx2).*slon(:,1),[3 1 2]);
                idx1 = sub2ind([nlat nlon nt], ilat+1 , ilons(:,1) , min(it+1,size(data,3)));
                idx2 = sub2ind([nlat nlon nt], ilat+1 , ilone(:,1) , min(it+1,size(data,3)));
                valad = permute(data(idx1).*(1-slon(:,1)) + data(idx2).*slon(:,1),[3 1 2]);

                % interpolate along lat
                % before
                % after
                valb = valbu(:).*(1-slat) + valbd(:).*slat;
                vala = valau(:).*(1-slat) + valad(:).*slat;

                % interpolate along time
                val = valb.*(1-st) + vala.*st;
            end
            
        end
        
        function val = linInterpRotLatLonTime(data, first_lat, dlat, first_lon, dlon, first_t, dt, lat, lon, t)
            % Interpolate values froma data on a geographical grid with multiple epoch
            % data structure:
            %        first dimension : dlat (+) south pole -> north pole
            %        second dimension : dlon (+) west -> east
            %        third dimension : dr (+) time usual direction
            %        NOTE: dlat, dlon,dt do not have to be positive
            %
            % Differently to linInterpLatLonTime it consider that the map rotates with the Earth motion
            % This function is created ad hoc for ionospheric interpolation
            %
            % INPUT:
            %      data - the data to be interpolate
            %      fist_lat - value of first lat value (max lat)
            %      dlat - px size lat
            %      first_lon - value of first lon value
            %      dlon - px size lon
            %      first_t - value of first time
            %      dt - px size time
            %      lat - lat at what we want to interpolate
            %      lon - lon at what we ant to interpolate
            %      gps_time - time at what we want to interpolate
            % NOTES 1 - all lat values should have same unit of measure
            %       2 - all lon values should have same unit of measure
            %       3 - all time values should have same unit of measure
            %       4 - the method will interpolate first in the dimesnion with less time
            % IMPORTANT : no double values at the borders should coexist: e.g. -180 180 or 0 360
            [nlat , nlon, nt] = size(data);
            n_in_lat = length(lat);
            n_in_lon = length(lon);
            n_in_t = length(t);
            assert(n_in_lat == n_in_lon);
            [ it, st, ilons, ilone, slon, ilat, slat] = Core_Utils.getIntIdxRot(data, first_lat, dlat, first_lon, dlon, first_t, dt, lat, lon, t);
            
            % lon first % NOTE: consider speed up in case only one time is present, unnecessary operations done
            % interpolate along lon 
            % before up
            % after up
            % before down
            % after down
            it = it*ones(size(ilat));            
            idx1 = sub2ind([nlat nlon nt], ilat, ilons(:,1), it);
            idx2 = sub2ind([nlat nlon nt], ilat, ilone(:,1), it);
            valbu = permute(data(idx1).*(1-slon(:,1)) + data(idx2).*slon(:,1),[3 1 2]);
            idx1 = sub2ind([nlat nlon nt], ilat   , ilons(:,2) , min(it+1,size(data,3)));
            idx2 = sub2ind([nlat nlon nt], ilat   , ilone(:,2) , min(it+1,size(data,3)));
            valau = permute(data(idx1).*(1-slon(:,2)) + data(idx2).*slon(:,2),[3 1 2]);
            idx1 = sub2ind([nlat nlon nt], ilat+1 , ilons(:,1) , it);
            idx2 = sub2ind([nlat nlon nt], ilat+1 , ilone(:,1) , it);
            valbd = permute(data(idx1).*(1-slon(:,1)) + data(idx2).*slon(:,1),[3 1 2]);
            idx1 = sub2ind([nlat nlon nt], ilat+1 , ilons(:,2) , min(it+1,size(data,3)));
            idx2 = sub2ind([nlat nlon nt], ilat+1 , ilone(:,2) , min(it+1,size(data,3)));
            valad = permute(data(idx1).*(1-slon(:,2)) + data(idx2).*slon(:,2),[3 1 2]);

            % interpolate along lat
            % before
            % after
            valb = valbu(:).*(1-slat) + valbd(:).*slat;
            vala = valau(:).*(1-slat) + valad(:).*slat;

            % interpolate along time
            val = valb.*(1-st) + vala.*st;            
        end
        
        function [dlat, dlon, data] = addSphBorder(dlat, dlon, data, border_size)
            % Replicate spherical data to allow interpolation
            %
            % INPUT
            %   dlat    latitude  [deg]  [  90 -90]
            %   dlon    longitude [deg]  [-180 180]
            %   data    data to replicate
            %
            % OUTPUT
            %   !!! Angles will be outside their range of validity !!!
            %   dlat    latitude  [deg]
            %   dlon    longitude [deg]
            %   data    data to replicate
            %
            % SYNTAX
            %   [dlat, dlon, data] = Core_Utils.addSphBorder(dlat, dlon, data, border_size)

            dlon = mod(dlon + 180, 360) - 180;
            idb = abs(dlat(:)) > 90 - border_size(1) & abs(dlat(:)) < 90;
            dlat = [dlat(:); (180 - abs(dlat(idb))) .* sign(dlat(idb))];
            dlon = [dlon(:); dlon(idb)];
            data = [data(:); data(idb)];

            idb = abs(dlon(:)) > 180 - border_size(end) & abs(dlat(:)) < 180;
            dlat = [dlat(:); dlat(idb)];
            dlon = [dlon(:); (360 - abs(dlon(idb))) .* -sign(dlon(idb))];
            data = [data(:); data(idb)];
        end
        
        function [data_grid, dlat_grid, dlon_grid] = sparseEarthGridder(dlat, dlon, data, grid_step)
            % Grid from spherical sparse data to the entire world map
            %
            % INPUT 
            %   dlat    latitude  [deg]  [  90 -90] [n x 1]
            %   dlon    longitude [deg]  [-180 180] [n x 1]
            %   data    data to interpolate         [n x 1]
            %
            % OUTPUT 
            %   data         Earth map 
            %   dlat_grid    latitude  [deg]  [  90 -90]
            %   dlon_grid    longitude [deg]  [-180 180]
            %
            % SYNTAX
            %   [data_grid, dlat_grid, dlon_grid] = Core_Utils.sparseEarthGridder(dlat, dlon, data, grid_step)


            border_size = 15;
            if nargin == 3 || isempty(grid_step)
                grid_step = [0.5 0.5];
            end

            [dlat, dlon, data] = Core_Utils.addSphBorder(dlat, dlon, data, border_size);

            fun = scatteredInterpolant(dlon(:), dlat(:), data(:));

            [dlat_grid, dlon_grid] = getGrid(grid_step);
            [dlon_grid, dlat_grid] = meshgrid(dlon_grid,dlat_grid);
            data_grid = reshape(fun(dlon_grid,dlat_grid), size(dlat_grid,1), size(dlat_grid,2));
        end

        function [val] = cubicSpline(t)
            % Compute matrix entry for cubic spline
            %
            % INPUT
            %   t -> 0 : 1
            %   order -> 1,3
            %
            % SYNTAX:
            %  Core_Utils.cubicSplic(t)
            val = zeros(numel(t),4);
            val(:,1) = (1 - t).^3/6;
            val(:,2) = ((2-t).^3 - 4*(1-t).^3)/6;
            val(:,3) = ((1+t).^3 - 4*(t).^3)/6;
            val(:,4) = (t).^3/6;
        end
        
        function [val] = linearSpline(t)
            % Compute matrix entry for linear spline
            %
            % INPUT
            %   t -> 0 : 1
            %   order -> 1,3
            %
            % SYNTAX:
            %  Core_Utils.cubicSplic(t)
            val = zeros(numel(t),2);
            val(:,1) = 1 -t;
            val(:,2) = t;
        end
        
        function [idx, val] = hemisphereCubicSpline(n_az, n_el, az, el)
            % give the index of the hemisphere spline idx
            % first the equator then the first parallel then the second
            % parallel
            %
            % SYNTAX:
            %  [idx] = hemisphereSpline(n_az,n_el,az,el)
            el_step = 90/n_el;
            idx_el = repmat(ceil(el / el_step),1,4);
            idx_el(:,2) = idx_el(:,2) + 1;
            idx_el(:,3) = idx_el(:,3) + 2;
            idx_el(:,4) = idx_el(:,4) + 3;
            az_step = 360/n_el;
            idx_az = repmat(ceil(az / az_step),1,4);
            idx_az(:,2) = idx_az(:,2) + 1;
            idx_az(:,3) = idx_az(:,3) + 2;
            idx_az(:,4) = idx_az(:,4) + 3;
            idx_az(idx_az > n_az) = idx_az(idx_az > n_az) - n_az;
            idx = idx_az + (idx_el -1).*n_el;
            idx = [idx idx+1 idx+2 idx+3];
            
            t_el = rem(el/el_step);
            t_az = rem(az/az_step);
            val_el = Core_Utils.cubicSpline(t_el);
            val_az = Core_Utils.cubicSpline(t_az);
            
            val = [val_az.*repmat(val_el(:,1),1,4) val_az.*repmat(val_el(:,2),1,4) val_az.*repmat(val_el(:,3),1,4) val_az.*repmat(val_el(:,4),1,4)];
        end
                
        function x = despline(x, spline_base)
            % despline signal x
            %
            % SYNTAX:
            %   x = despline(x,<spline_base>)
            if nargin <2
                spline_base = round(size(x,1)/7);
            end
            for i = 1: size(x,2)
                x(:,2) = x(:,2) - splinerMat(1:length(x(:,2)),x(:,2),spline_base);
            end
        end
        
        
        function [az_grid, el_grid] = getPolarGrid(step_az, step_el)
            % Get a regularly spaced lnots coordinates for a polar grid
            %
            % INPUT
            %   step_az     step of the grid in azimuth   [deg]
            %   step_el     step of the grid in elevation [deg]
            %
            % OUTPUT 
            %   el_grid     grid centers in azimuth   [deg]
            %   el_grid     grid centers in elevation [deg]
            %
            % SYNTAX
            %   [az_grid, el_grid] = Core_Utils.getPolarGrid(step_az, step_el)
            [el_grid, az_grid] = getGrid([step_el step_az], 0, 90, -180, 180);
        end
            
        function id_ok = polarCleaner(az, el, data, step_deg, n_sigma)
            % Remove observations above n_sigma sigma 
            %
            % INPUT 
            %   az          azimuth [rad]
            %   el          elevation [rad]
            %  data
            %  step_deg     horizontal array with max two lines [2x2]
            %               step_deg(1,:) size of cells - flag everything below 3 sigma 
            %               step_deg(2,:) size of cells - unflag everything below 3 sigma 
            %
            % SYNTAX
            %   id_ok = Core_Utils.polarCleaner(az, el, data, step_deg)
            
            % az -180 : 180
            % el 0 : 90
            %figure; polarScatter(az, pi/2-el, 5, data); colormap(gat);
            if nargin < 4
                step_deg = [360, 1; 3, 3];
            end
            if nargin < 5
                n_sigma = 3;
            end
            
            
            az_grid = ((-180 + (step_deg(1, 1) / 2)) : step_deg(1, 1) : (180 - step_deg(1, 1) / 2)) .* (pi/180);
            el_grid = flipud(((step_deg(1, end) / 2) : step_deg(1, end) : 90 - (step_deg(1, end) / 2))' .* (pi/180));
            n_az = numel(az_grid);
            n_el = numel(el_grid);
            
            % Find map indexes
            col = max(1, min(floor((az + pi) / (step_deg(1, 1) / 180 * pi) ) + 1, length(az_grid)));
            row = max(1, min(floor((pi/2 - el) / (step_deg(1, end) / 180 * pi)) + 1, length(el_grid)));
            
            uid = row + (col-1) * numel(el_grid);
            id_ok = true(size(data));
            for b = unique(uid)'
                id_set = find(uid == b);
                dset = data(id_set);
                sigma = std(dset);
                id_ok(id_set(abs(dset) > n_sigma * sigma)) = false;
            end
            %figure; polarScatter(az(id_ok == 0), pi/2-el(id_ok == 0), 5, data(id_ok == 0)); colormap(gat);
            
            if size(step_deg, 1) == 2
                step_deg = step_deg(2,:);
                
                az_grid = ((-180 + (step_deg(1, 1) / 2)) : step_deg(1, 1) : (180 - step_deg(1, 1) / 2)) .* (pi/180);
                el_grid = flipud(((step_deg(1, end) / 2) : step_deg(1, end) : 90 - (step_deg(1, end) / 2))' .* (pi/180));
                n_az = numel(az_grid);
                n_el = numel(el_grid);
                
                % Find map indexes
                col = max(1, min(floor((az + pi) / (step_deg(1) / 180 * pi) ) + 1, length(az_grid)));
                row = max(1, min(floor((pi/2 - el) / (step_deg(end) / 180 * pi)) + 1, length(el_grid)));
                
                uid = row + (col-1) * numel(el_grid);
                for b = unique(uid)'
                    id_set = find(uid == b);
                    dset = data(id_set) - median(serialize(data(id_set)));
                    %data(id_set) = median(serialize(data(id_set)));
                    sigma = min(0.01, std(dset));
                    id_ok(id_set(abs(dset) < n_sigma * sigma)) = true;
                end
                
                %figure; polarScatter(az(id_ok), pi/2-el(id_ok), 5, data(id_ok)); colormap(gat);
                %figure; polarScatter(az(id_ok == 0), pi/2-el(id_ok == 0), 5, data(id_ok == 0)); colormap(gat);
            end
            
        end
        
               
        function [data_map, n_map, az_grid_out, el_grid_out] = hemiGridder(az, el, data, step_deg, step_deg_out, flag_congurent_cells, n_min)
            % Grid points on a regularly gridded semi sphere
            %
            % hemiGridder grids points on a regularly gridded semi-sphere.
            % This function is used for spatially gridding data points defined by azimuth
            % and elevation onto a semi-spherical surface. It supports both uniform and
            % non-uniform gridding, and can handle sparse data distributions.
            %
            % INPUTS:
            %   az                      - Azimuth values [degrees]. Range: -180 to 180.
            %   el                      - Elevation values [degrees]. Range: 0 to 90.
            %   data                    - Data values corresponding to each az, el pair.
            %   step_deg                - Gridding step size in degrees for uniform gridding.
            %   step_deg_out            - Output step size in degrees for re-gridding.
            %   flag_congurent_cells    - Boolean flag to enable non-uniform gridding based on
            %                             congruent cells (true) or use uniform gridding (false).
            %   n_min                   - Minimum number of data points required in a cell.
            %                             Cells with fewer points will be set to zero.
            %
            % OUTPUTS:
            %   data_map                - 2D grid of data values on the semi-sphere.
            %   n_map                   - 2D grid representing the count of data points in each cell.
            %   az_grid_out             - Grid azimuth values post re-gridding.
            %   el_grid_out             - Grid elevation values post re-gridding.
            %
            % SYNTAX:
            %   [data_map, n_map, az_grid_out, el_grid_out] =
            %       hemiGridder(az, el, data, step_deg, step_deg_out, flag_congurent_cells, n_min)
            %
            % DETAILS:
            %   The function first defines the grid based on the input parameters, then
            %   maps the data onto this grid. If the flag_congurent_cells is set, the azimuth
            %   gridding varies with elevation to maintain roughly equal area cells across the
            %   semi-sphere. The function also supports re-gridding the data to a different
            %   resolution specified by step_deg_out. Post-processing options include filtering
            %   out cells with insufficient data points.
            %
            % EXAMPLE:
            %   [data_map, n_map, az_grid, el_grid] =
            %       hemiGridder(az, el, data, 5, 2, true, 3);
            %   This will grid the data on a semi-sphere with 5-degree gridding, re-grid with 2-degree
            %   resolution, use non-uniform gridding, and filter out cells with fewer than 3 data points.
            %
            el_grid = flipud(((step_deg(end) / 2) : step_deg(end) : 90 - (step_deg(end) / 2))' .* (pi/180));
            flag_congurent_cells = nargin >= 6 && ~isempty(flag_congurent_cells) && flag_congurent_cells;
            if nargin < 7
                n_min = 0; % by default grid all the data
            end            
            if flag_congurent_cells
                step_az = 360 ./ round((360 / step_deg(1)) * cos(el_grid));
                az_grid = {};
                for i = 1 : numel(step_az)
                   az_grid{i} = ((-180 + (step_az(i) / 2)) : step_az(i) : (180 - step_az(i) / 2)) .* (pi/180);
                   n_az(i) = numel(az_grid{i});
                end                
            else
                az_grid = ((-180 + (step_deg(1) / 2)) : step_deg(1) : (180 - step_deg(1) / 2)) .* (pi/180);
                n_az = numel(az_grid);
            end
            n_el = numel(el_grid);
            
            % Find map indexes
            row = max(1, min(floor((pi/2 - el) / (step_deg(end) / 180 * pi)) + 1, length(el_grid)));
            if flag_congurent_cells
                col = max(1, min(floor((az + pi) ./ (step_az(row) / 180 * pi) ) + 1, n_az(i)));
            else
                col = max(1, min(floor((az + pi) / (step_deg(1) / 180 * pi) ) + 1, length(az_grid)));
            end
            
            % init maps
            [n_map, data_map] = deal(zeros(n_el, max(n_az)));
                        
            % fill maps
            for i = 1 : numel(data)
                n_map(row(i), col(i)) = n_map(row(i), col(i)) + 1;
                data_map(row(i), col(i)) = data_map(row(i), col(i)) + data(i);
            end
            data_map(n_map > 0) = data_map(n_map > 0) ./ (n_map(n_map > 0));
            
            % zeros cells with a minimum number of data < n_min
            n_map = nan2zero(n_map);
            data_map(n_map <= n_min) = 0;

            flag_debug = false;
            if flag_congurent_cells && flag_debug
                data_congruent = {};
                for i = 1 : numel(el_grid)
                    data_congruent{i} = data_map(i, 1 : numel(az_grid{i}));
                end
                Core_Utils.plotSphPatchGrid(el_grid, az_grid, data_congruent);
            end
            
            % distort map
            if (nargin >= 5) && ~isempty(step_deg_out)
                data_map_in = data_map;
                decl_n = ((pi/2 - el_grid)/(pi/2));
                
                % Define output grid
                az_grid_out = ((-180 + (step_deg(1) / 2)) : step_deg(1) : (180 - step_deg(1) / 2)) .* (pi/180);
                el_grid_out = el_grid;
                n_az_out = numel(az_grid_out);
                n_el_out = numel(el_grid_out);
                data_map = zeros(n_el_out, n_az_out);
                
                if flag_congurent_cells
                    % Interpolate elevation by elevation
                    [az, el] = deal(zeros(n_el, max(n_az)));
                    [az_mg, el_mg] = meshgrid(az_grid_out, el_grid_out);
                    for i = 1 : numel(el_grid)
                        az_tmp = [az_grid{i} nan(1, max(n_az) - n_az(i))];
                        az(i, :) = az_tmp;
                        el(i, :) = el_grid(i);
                        
                        if sum(n_map(i, :) > 0) < 2
                            data_map(i, :) = data_map_in(i,1);
                        else
                            az_tmp = az_tmp(n_map(i, :) > 0)';
                            az_tmp = [az_tmp-2*pi; az_tmp; az_tmp+2*pi];
                            data_tmp = data_map_in(i, n_map(i, :) > 0)';
                            data_tmp = [data_tmp; data_tmp; data_tmp];
                            data_map(i, :) = interp1(az_tmp, data_tmp, az_mg(i, :)', 'linear');
                            az_tmp = [az_grid{i}'-2*pi; az_grid{i}'; az_grid{i}'+2*pi];
                            data_tmp = [n_map(i, 1 : n_az(i))'; n_map(i, 1 : n_az(i))'; n_map(i, 1 : n_az(i))'];
                            n_map(i, :) = interp1(az_tmp, data_tmp, az_mg(i, :)', 'nearest');
                        end
                    end
                    data_map(n_map == 0) = 0;
                    n_map(n_map == 0) = 0.1; % this is to cheat the next scatteredInterpolant
                    data_map_in = data_map;
                    
                    % get polar coordinates
                    x = sin(az_grid_out) .* decl_n;
                    y = cos(az_grid_out) .* decl_n;
                    funGridder = scatteredInterpolant(x(n_map > 0), y(n_map > 0), data_map_in(n_map > 0), 'linear' );
                else
                    % get polar coordinates
                    x = sin(az_grid) .* decl_n;
                    y = cos(az_grid) .* decl_n;
                    
                    data_map_in(n_map == 0) = 0;
                    n_map(n_map == 0) = 0.1; % this is to cheat the next scatteredInterpolant
                    
                    funGridder = scatteredInterpolant(x(n_map > 0), y(n_map > 0), data_map_in(n_map > 0), 'linear' );
                    x = linspace(-1, 1, 180/min(step_deg_out));
                    y = x;
                    [x_mg, y_mg] = meshgrid(x, y);
                    polar_data = nan(numel(x), numel(y));
                    id_ok = hypot(x_mg, y_mg) < 1;
                    polar_data(id_ok) = funGridder(x_mg(id_ok), y_mg(id_ok));
                    
                    % Prepare polar gridder
                    funGridder = scatteredInterpolant(x_mg(id_ok), y_mg(id_ok), polar_data(id_ok), 'linear');
                end
                
                
                % Define output grid
                az_grid_out = ((-180 + (step_deg_out(1) / 2)) : step_deg_out(1) : (180 - step_deg_out(1) / 2)) .* (pi/180);
                el_grid_out = flipud(((step_deg_out(end) / 2) : step_deg_out(end) : 90 - (step_deg_out(end) / 2))' .* (pi/180));
                n_az_out = numel(az_grid_out);
                n_el_out = numel(el_grid_out);
                data_map = zeros(n_el_out, n_az_out);
                
                % Get polar coordinates
                decl_n = ((pi/2 - el_grid_out)/(pi/2));
                x = sin(az_grid_out) .* decl_n;
                y = cos(az_grid_out) .* decl_n;
                
                data_map(:) = funGridder(x(:),y(:));
            else
                el_grid_out = el_grid;
                az_grid_out = az_grid;
            end
        end
        
        function [data_map, n_map, az_grid_out, el_grid_out] = hemiGridderStd(az, el, data, step_deg, step_deg_out, flag_congurent_cells)
            % Get std of points on a semi sphere
            % To be fixed for non congruent grids
            %
            % INPUT 
            %   az      azimuth
            %   el      elevation
            %
            % SYNTAX
            %   [data_map, n_map, az_grid, el_grid] = Core_Utils.hemiGridderStd(az, el, data, step_deg, step_deg_out, flag_congurent_cells)
            
            % Define grid
            % az -180 : 180
            % el 0 : 90
            el_grid = flipud(((step_deg(end) / 2) : step_deg(end) : 90 - (step_deg(end) / 2))' .* (pi/180));
            flag_congurent_cells = nargin >= 6 && ~isempty(flag_congurent_cells) && flag_congurent_cells;
            if flag_congurent_cells
                step_az = 360 ./ round((360 / step_deg(1)) * cos(el_grid));
                az_grid = {};
                for i = 1 : numel(step_az)
                   az_grid{i} = ((-180 + (step_az(i) / 2)) : step_az(i) : (180 - step_az(i) / 2)) .* (pi/180);
                   n_az(i) = numel(az_grid{i});
                end                
            else
                az_grid = ((-180 + (step_deg(1) / 2)) : step_deg(1) : (180 - step_deg(1) / 2)) .* (pi/180);
                n_az = numel(az_grid);
            end
            n_el = numel(el_grid);
            
            % Find map indexes
            row = max(1, min(floor((pi/2 - el) / (step_deg(end) / 180 * pi)) + 1, length(el_grid)));
            if flag_congurent_cells
                col = max(1, min(floor((az + pi) ./ (step_az(row) / 180 * pi) ) + 1, n_az(i)));
            else
                col = max(1, min(floor((az + pi) / (step_deg(1) / 180 * pi) ) + 1, length(az_grid)));
            end
            
            % init maps
            [n_map, data_map] = deal(zeros(n_el, max(n_az)));
                        
            % fill maps
            r_lim = minMax(row);
            if any(r_lim)
                for r = r_lim(1):r_lim(2)
                    id_r = row == r;
                    for c = unique(col(id_r))'
                        id = id_r & col == c; % index of the cell
                        n_map(r,c) = sum(id);
                        data_map(r,c) = std(data(id));
                    end
                end
            end
            
            % max cells with a minimum number of data < n_min
            n_map = nan2zero(n_map);
            data_map(n_map < 3) = min(max(data_map(:)), 3*std(data(:),'omitnan'));

            flag_debug = false;
            if flag_congurent_cells && flag_debug
                data_congruent = {};
                for i = 1 : numel(el_grid)
                    data_congruent{i} = data_map(i, 1 : numel(az_grid{i}));
                end
                Core_Utils.plotSphPatchGrid(el_grid, az_grid, data_congruent);
            end
            
            % distort map
            if (nargin >= 5) && ~isempty(step_deg_out)
                data_map_in = data_map;
                decl_n = ((pi/2 - el_grid)/(pi/2));
                
                % Define output grid
                az_grid_out = ((-180 + (step_deg(1) / 2)) : step_deg(1) : (180 - step_deg(1) / 2)) .* (pi/180);
                el_grid_out = el_grid;
                n_az_out = numel(az_grid_out);
                n_el_out = numel(el_grid_out);
                data_map = zeros(n_el_out, n_az_out);
                
                if flag_congurent_cells
                    % Interpolate elevation by elevation
                    [az, el] = deal(zeros(n_el, max(n_az)));
                    [az_mg, el_mg] = meshgrid(az_grid_out, el_grid_out);
                    for i = 1 : numel(el_grid)
                        az_tmp = [az_grid{i} nan(1, max(n_az) - n_az(i))];
                        az(i, :) = az_tmp;
                        el(i, :) = el_grid(i);
                        
                        if sum(n_map(i, :) > 0) < 2
                            data_map(i, :) = data_map_in(i,1);
                        else
                            az_tmp = az_tmp(n_map(i, :) > 0)';
                            az_tmp = [az_tmp-2*pi; az_tmp; az_tmp+2*pi];
                            data_tmp = data_map_in(i, n_map(i, :) > 0)';
                            data_tmp = [data_tmp; data_tmp; data_tmp];
                            data_map(i, :) = interp1(az_tmp, data_tmp, az_mg(i, :)', 'linear');
                            az_tmp = [az_grid{i}'-2*pi; az_grid{i}'; az_grid{i}'+2*pi];
                            data_tmp = [n_map(i, 1 : n_az(i))'; n_map(i, 1 : n_az(i))'; n_map(i, 1 : n_az(i))'];
                            n_map(i, :) = interp1(az_tmp, data_tmp, az_mg(i, :)', 'nearest');
                        end
                    end
                    %data_map(n_map == 0) = 0;
                    %n_map(n_map == 0) = 0.1; % this is to cheat the next scatteredInterpolant
                    data_map_in = data_map;
                    
                    % get polar coordinates
                    x = sin(az_grid_out) .* decl_n;
                    y = cos(az_grid_out) .* decl_n;
                    funGridder = scatteredInterpolant(x(n_map > 0), y(n_map > 0), data_map_in(n_map > 0), 'linear' );
                    funGridderNMap = scatteredInterpolant(x(:), y(:), n_map(:), 'nearest' );
                else
                    % get polar coordinates
                    x = sin(az_grid) .* decl_n;
                    y = cos(az_grid) .* decl_n;
                    
                    data_map_in(n_map == 0) = 0;
                    n_map(n_map == 0) = 0.1; % this is to cheat the next scatteredInterpolant
                    
                    funGridder = scatteredInterpolant(x(n_map > 0), y(n_map > 0), data_map_in(n_map > 0), 'linear' );
                    x = linspace(-1, 1, 180/min(step_deg_out));
                    y = x;
                    [x_mg, y_mg] = meshgrid(x, y);
                    polar_data = nan(numel(x), numel(y));
                    id_ok = hypot(x_mg, y_mg) < 1;
                    polar_data(id_ok) = funGridder(x_mg(id_ok), y_mg(id_ok));
                    
                    % Prepare polar gridder
                    funGridder = scatteredInterpolant(x_mg(id_ok), y_mg(id_ok), polar_data(id_ok), 'linear');
                end
                
                
                % Define output grid
                az_grid_out = ((-180 + (step_deg_out(1) / 2)) : step_deg_out(1) : (180 - step_deg_out(1) / 2)) .* (pi/180);
                el_grid_out = flipud(((step_deg_out(end) / 2) : step_deg_out(end) : 90 - (step_deg_out(end) / 2))' .* (pi/180));
                n_az_out = numel(az_grid_out);
                n_el_out = numel(el_grid_out);
                data_map = zeros(n_el_out, n_az_out);
                n_map_out = zeros(n_el_out, n_az_out);
                
                % Get polar coordinates
                decl_n = ((pi/2 - el_grid_out)/(pi/2));
                x = sin(az_grid_out) .* decl_n;
                y = cos(az_grid_out) .* decl_n;
                
                data_map(:) = funGridder(x(:),y(:));
                n_map_out(:) = funGridderNMap(x(:),y(:));
                data_lim = minMax(data_map(n_map_out > 3));
                data_map(data_map < data_lim(1)) = data_lim(1);
                data_map(data_map > data_lim(2)) = data_lim(2);
                for a = 1 : size(data_map,2)
                    data_map(find(n_map_out(:,a) > 1, 1, 'last') : end, a) = data_lim(2);
                end
                n_map = n_map_out;
            else
                el_grid_out = el_grid;
                az_grid_out = az_grid;
            end
        end
        
        function [data, row, col] = hgrid2scatter(az, el, hmap, flag_congurent_cells, method)
            % Grid points on a regularly gridded semi sphere
            %
            % INPUT 
            %   az      azimuth list   [n x 1]
            %   el      elevation list [n x 1]
            %   hmap    gridded hemispherical map
            %   
            %
            % SYNTAX
            %   [data, row, col] = Core_Utils.hgrid2scatter(az, el, hmap)
            
            % Define grid
            % az -180 : 180
            % el 0 : 90
            
            step_deg = fliplr([90 360] ./ size(hmap));
            el_grid = flipud(((step_deg(end) / 2) : step_deg(end) : 90 - (step_deg(end) / 2))' .* (pi/180));
            flag_congurent_cells = nargin >= 4 && ~isempty(flag_congurent_cells) && flag_congurent_cells;
            
            if nargin < 5 || isempty(method)
                method = 'spline';
            end            
            if flag_congurent_cells
                method = 'nearest';
            end
            if flag_congurent_cells
                step_az = 360 ./ round((360 / step_deg(1)) * cos(el_grid));
                az_grid = {};
                for i = 1 : numel(step_az)
                   az_grid{i} = ((-180 + (step_az(i) / 2)) : step_az(i) : (180 - step_az(i) / 2)) .* (pi/180);
                   n_az(i) = numel(az_grid{i});
                end                
            else
                az_grid = ((-180 + (step_deg(1) / 2)) : step_deg(1) : (180 - step_deg(1) / 2)) .* (pi/180);
                n_az = numel(az_grid);
            end
            n_el = numel(el_grid);
            
            if strcmp(method, 'nearest')
                % Find map indexes
                row = max(1, min(floor((pi/2 - el) / (step_deg(end) / 180 * pi)) + 1, length(el_grid)));
                if flag_congurent_cells
                    col = max(1, min(floor((az + pi) ./ (step_az(row) / 180 * pi) ) + 1, n_az(i)));
                else
                    col = max(1, min(floor((az + pi) / (step_deg(1) / 180 * pi) ) + 1, length(az_grid)));
                end
                
                % init data
                data = nan(size(az));
                
                % fill maps
                for i = 1 : numel(data)
                    data(i) = hmap(row(i), col(i));
                end
            else
                % Add repeating padding
                n_col = 3;
                az_grid = [(az_grid(end-n_col+1:end) - 2*pi) az_grid (az_grid(1:n_col) + 2*pi)];
                hmap = [hmap(:, end-n_col+1:end) hmap hmap(:, 1:n_col)];
                % Interp
                [a, e] = ndgrid(az_grid', flipud(el_grid));
                finterp = griddedInterpolant(a,e,flipud(hmap)', method);
                data = finterp(az, el);
            end
        end
                
                
        function [n] = getAllNormCoeffZernike(l_max, m_max)
            % Generate all the Zernike normalization coefficient
            %
            % SINTAX
            %   [n] = getAllNormCoeffZernike(l_max, m_max)
            n_par = l_max * (l_max + 3) / 2 + 1;
            l = zeros(n_par, 1);
            m = zeros(n_par, 1);
            i = 0;
            for degree = 0 : l_max
                i = i(end) + (1 : degree + 1);
                l(i) = degree;
                m(i) = -degree : 2 : degree;
            end
            l(abs(m) > m_max) = [];
            m(abs(m) > m_max) = [];
            n = sqrt((1+(m~=0)).*(l+1)/pi);
        end
        
        function [l , m] = getAllZdegree(l_max, m_max)
            % Generate all the Zernike degree combinations
            %
            % INPUT
            %   l_max       max degrees
            %   m_max       max orders
            
            n_par = l_max * (l_max + 3) / 2 + 1;
            
            l = zeros(n_par, 1);
            m = zeros(n_par, 1);
            i = 0;
            for degree = 0 : l_max
                i = i(end) + (1 : degree + 1);
                l(i) = degree;
                m(i) = -degree : 2 : degree;
            end
            
            l(abs(m) > m_max) = [];
            m(abs(m) > m_max) = [];
        end
        
        function [z, l, m] = getAllZernikeNorm(l_max, m_max, az, r)
            % Generate all the Zernike parameters combinations
            %
            % INPUT
            %   l       list of degrees
            %   m       list of orders
            %   az      list of azimuth angles   [rad]
            %   r       disc radius              [0 1]
            %
            % SINTAX
            %   [z, l, m] = getAllZernikeNorm(l_max, az, el)
            %   [z, l, m] = getAllZernikeNorm(l_max, m_max, az, el)
            if nargin == 3
                r = az;
                az = m_max;
                m_max = l_max;
            end
            
            n_par = l_max * (l_max + 3) / 2 + 1;
            
            l = zeros(n_par, 1);
            m = zeros(n_par, 1);
            i = 0;
            for degree = 0 : l_max
                i = i(end) + (1 : degree + 1);
                l(i) = degree;
                m(i) = -degree : 2 : degree;
            end
            
            l(abs(m) > m_max) = [];
            m(abs(m) > m_max) = [];
            
            z = Core_Utils.getZernikeNorm(l, m, az, r);
        end
        
        function z = getZernikeNorm(l, m, az, r)
            % Get Zernike values for the polynomials
            %
            % INPUT 
            %   l       list of degrees
            %   m       list of orders
            %   az      list of azimuth angles   [rad]
            %   r       disc radius              [0 1]
            %
            % SINTAX
            %   z = getZernike(l, m, az, el)
            
            z = zernfun(l, m, r(:), az(:), 'norm');
        end
                        
        function [S] = reorthZernikeMask(lat, lon, el_thrsh, n_sample)
            % get an orthogonal basis of zernike function that maximize the
            % power in the area of interest
            %
            % SYNTAX:
            % [S] = Core_Utils.reorthZernikeMask(lat,lon,el_thrsh)
            if nargin < 4
                n_sample = 500000;
            end
            %             N = zeros(741,741);
            %             for j = 1 :10
            xy = (rand(n_sample*2,2)-0.5)*2;
            len= sqrt(xy(:,1).^2 + xy(:,2).^2);
            xy(len > 1,:) = [];
            len(len > 1) = [];
            r = len;
            theta = atan2(xy(:,2),xy(:,1));
            el = pi/2 - r*pi/2;
            cc = Constellation_Collector();
            xy(el < (el_thrsh-(2/180*pi)),:) = []; %remove under treshold
            [mask_north, mask_sud] = cc.getGPS.getPolarMask(lat,lon,5);
            
            xy_mask_north = [sin(mask_north(:,1)).*(1 - mask_north(:,2)/(pi/2)) cos(mask_north(:,1)).*(1 - mask_north(:,2)/(pi/2)) ];
            xy_mask_north = [xy_mask_north; xy_mask_north(1,:)];
            idx_inside = true(size(xy,1),1);
            for i = 1 : (size(xy_mask_north,1)-1)
                idx_inside = idx_inside & ((xy_mask_north(i+1,1) - xy_mask_north(i,1)) * (xy(:,2) - xy_mask_north(i,2)) - (xy_mask_north(i+1,2) - xy_mask_north(i,2)) * (xy(:,1) - xy_mask_north(i,1))) > 0;
            end
            xy(idx_inside,:) = [];
            
            len= sqrt(xy(:,1).^2 + xy(:,2).^2);
            r = len;
            theta = atan2(xy(:,2),xy(:,1));
            el = pi/2 - r*pi/2;
            
            [z1] = Zernike.getAll(10, 10, theta, r);
            %z1 = [zernfun(1,1,r,theta,'norm')      zernfun(1,-1,r,theta,'norm') zernfun(5,5,r,theta,'norm')];%           2/sqrt(pi)
            %z1 = [r .* sin(theta)      r .* cos(theta)];
            %       1   -1    r * sin(theta)                 2/sqrt(pi)
            
            %             N = N + z1'*z1;
            %             j
            %             end
            N = z1'*z1;
            [U,D,V] = svd(N);
            d = diag(D);
            % get the flexum
            exclude_fisrt = round(length(d)/5);
            [power_smooth] = splinerMat(exclude_fisrt:length(d),d(exclude_fisrt:end),20);
            [~,flex] = max(Core_Utils.diffAndPred(power_smooth,3));
            
            stop = exclude_fisrt + flex;
            
            
            S = U(:,1:stop);
        end
        
        
        function flag_is_hold = isHold()
            % Return if there is an open figure on hold
            %
            % SYNTAX
            %   Core_Utils.isHold()
            
            flag_is_hold = ~isempty(findobj('Type', 'figure')) && ishold;
        end
        
        function fh = showZernike3(l, m, z_par, el_min)
            % Show 3D plot of Zernike polynomials 
            %
            % SINTAX
            %   fh = showZernike3(l, m, z_par)
            
            % [x, y] = pol2cart(theta, r_synt);
            if nargin == 4
                r_max = 1 - (2 * el_min / pi);
                [theta, r_prj] = meshgrid(linspace(0, 2*pi, 361), linspace(0, r_max, 101));
            else
                [theta, r_prj] = meshgrid(linspace(0, 2*pi, 361), linspace(0, 1, 101));
            end
            % r_prj is the correct radius for my polar projection 
            r_zern = r_prj;
    
            z = nan(size(theta));
            z(:) = zernfun(l, m, r_zern(:), theta(:)) * z_par;
            %z = r_zern.*cos(theta);
            
            
            fh = figure();
            title('Zernike expansion')
            %polarplot3d(z, 'PlotType','surfn');
            polarplot3d(z,'PlotType','surfn','PolarGrid',{4 24}, 'PolarDirection', 'cw', 'TickSpacing',5,...
                   'RadLabels',4, 'RadLabelLocation',{180 'max'}, 'RadLabelColor','black', 'AxisLocation', 'mean');
            ar = get(gca,'DataAspectRatio');
            set(gca, 'DataAspectRatio', [1 1 ar(3)]);
            
            colormap(jet);
            material([ 0.4 0.9 0.55])
            l1 = light('position',[200 -300 400], 'color', [0.6 0.6 0.6]);
            l2 = light('position',[-600 600 900], 'color', [0.6 0.6 0.6]);
            l3 = light('position',[0 2 100], 'color', [0.6 0.6 0.6]);
            view(35, 45);
            Core_UI.beautifyFig(fh, 'dark');            
        end
        
        function fh = showZernike3StylePCV(l, m, z_par, el_min, limits)
            % Show 3D plot of Zernike polynomials 
            %
            % INPUT
            %   l       zernike degree [n x 1]
            %   m       zernike degree [n x 1]
            %   z_par   zernike degree [n x 1]
            %   el_min  cut-off angle  [rad]
            %   limits  min max (saturation) of the z map
            %
            % SINTAX
            %   fh = showZernike3(l, m, z_par, <el_min>, <limits>)
            
            % [x, y] = pol2cart(theta, r_synt);
            if nargin >= 4 && ~isempty(el_min)
                r_max = 1 - (2 * el_min / pi);
                [theta, r_prj] = meshgrid(linspace(0, 2*pi, 361), linspace(0, r_max, 101));
            else
                [theta, r_prj] = meshgrid(linspace(0, 2*pi, 361), linspace(0, 1, 101));
            end
            % r_prj is the correct radius for my polar projection 
            r_zern = r_prj;
    
            z = nan(size(theta));
            z(:) = zernfun(l, m, r_zern(:), theta(:)) * z_par;
            
            fh = figure();
            title('Zernike expansion')
            
            if nargin == 5 && ~isempty(limits)
                z = min(limits(2), max(z, limits(1)));
            end
            polarplot3d(z, 'PlotType', 'surf', 'RadialRange',[0 90] / 180 * pi, ...
                    'AxisLocation', 0, 'InterpMethod', 'cubic', ...
                    'PlotType', 'surfn', 'tickspacing', 15, ...
                    'GridColor', [0.7 0.7 0.7]);
                        
            colormap(flipud(Cmap.get('PuOr', 256)));
            
            axprop = {'DataAspectRatio',[1 1 8],'View', [-12 38], ...
                'Xlim', 1.5 * [-90 90] / 180 * pi, 'Ylim', 1.5 * [-90 90] / 180 * pi, ...
                % 'XTick', [], 'YTick', [], 'Color', 'none', 'XColor', 'none', 'YColor', 'none'
                };
            ax = gca;
            set(ax, axprop{:});
            
            if nargin >= 5 && ~isempty(limits)
                caxis(limits);
                zlim(limits);
            end

            
            %Core_UI.beautifyFig(fh, 'dark');
        end
        
        function fh = polarZerMap(l_max, m_max, az, r, data)
            % Take scattered observation and plot a polar interpolation
            %
            % INPUT
            %   l_max   maximum degree
            %   m_max   maximum order
            %   az      list of azimuth angles   [rad]
            %   r       disc radius              [0 1]
            %   data    list of data to be analized            
            %
            % SYNTAX
            %   fh = polarZerMap(l_max, m_max, az, r, data)
            [z_par, l, m] = Zernike.analysisAll(l_max, m_max, az, r, data, 1);
            fh = Zernike.showZernike(l, m, z_par);
        end
        
        function fh = polarZerMapDual(l_max, m_max, az, r, data)
            % Take scattered observation and plot:
            %  - a scattered polar
            %  - a polar interpolation
            %
            % INPUT
            %   l_max   maximum degree
            %   m_max   maximum order
            %   az      list of azimuth angles   [rad]
            %   r       disc radius              [0 1]
            %   data    list of data to be analized            
            %
            % SYNTAX
            %   fh = polarZerMapDual(l_max, m_max, az, el, data)
            fh = figure('Visible', 'off'); subplot(1,2,1); hold on;
            Core_Utils.polarZerMap(l_max, m_max, az, r, data);
            fh.Visible = false;
            cax = caxis();
            subplot(1,2,2); hold on; 
            polarScatter(az, (r * pi/2), 50, data, 'filled'); colorbar(); colormap(jet);
            cax_s = caxis();
            %cax = [min(cax(1), cax_s(1)) max(cax(2), cax_s(2))];
            cax = cax_s;
            caxis(cax);
            subplot(1,2,1);
            caxis(cax);
            fh.Visible = true;
            Core_UI.addExportMenu(fh); Core_UI.addBeautifyMenu(fh); Core_UI.beautifyFig(fh, 'dark');
        end
        
        function fh = polarZerMapQuad(l_max, m_max, az, r, data)
            % Take scattered observation and plot:
            %  - a scattered polar
            %  - a polar interpolation
            %  - a scattered vs el
            %  - an interpolation vs el
            %
            % INPUT
            %   l_max   maximum degree
            %   m_max   maximum order
            %   az      list of azimuth angles   [rad]
            %   r       disc radius              [0 1]
            %   data    list of data to be analized            
            %
            % SYNTAX
            %   fh = polarZerMapQuad(l_max, m_max, az, el, data)
            if ~isempty(findobj('Type', 'figure')) && ~ishold
                fh = gcf;
            else
                fh = figure('Visible', 'on'); 
            end
            
            subplot(2,2,1); hold on;
            
            Core_Utils.polarZerMap(l_max, m_max, az, r, data);
            cax = caxis();
            
            subplot(2,2,2); hold on; 
            polarScatter(az, r * pi/2, 50, data, 'filled'); colorbar(); colormap(jet);
            cax_s = caxis();
            %cax = [min(cax(1), cax_s(1)) max(cax(2), cax_s(2))];
            cax = cax_s;
            caxis(cax);
            subplot(2,2,1);
            caxis(cax);
            
            subplot(2,2,3); 
            data_smooth = Zernike.filter(l_max, m_max, az, r, data); hold on;       
            plot((r * pi/2) / pi * 180, data_smooth, '.', 'Color', Core_UI.getColor(1));
            
            el_grid = 0 : 90;
            plot(el_grid, Core_Utils.interp1LS((r * pi/2) / pi * 180, data_smooth, 3, el_grid), '--', 'LineWidth', 2, 'Color', Core_UI.getColor(3));
            
            ylim(cax);
            xlim([0 90]);
            grid minor;

            subplot(2,2,4); hold on; 
            plot((r * pi/2) / pi * 180, data, '.', 'Color', Core_UI.getColor(2)); hold on;
            plot(el_grid, Core_Utils.interp1LS((r * pi/2) / pi * 180, data, 3, el_grid), '--', 'LineWidth', 2 , 'Color', Core_UI.getColor(5));
            ylim(cax);
            xlim([0 90]);
            grid minor;
            
            fh.Visible = true;
            Core_UI.addExportMenu(fh); Core_UI.addBeautifyMenu(fh); Core_UI.beautifyFig(fh, 'dark');
        end
        
                
        function sphTest()
            %%  along track analysis            
            l_min = 0;
            l_max = 30;
            id = 1:20:length(data);
            [az el] = meshgrid(lambdaGrid, phiGrid);
            [ N_subSet, TN_subSet ] = analisiPointC_defineSYS (data(:), 1+0*data(:), el(:)./180*pi, az(:)./180*pi, el(:)*0+1, l_min, l_max, l_min, l_max, 1, 1, 0);
            tic;
            idM = 1:(l_max+1)^4;
            idM = reshape(idM,(l_max+1)^2,(l_max+1)^2);
            idM = tril(idM); idM = idM(idM>0);
            N = zeros((l_max+1)^2,(l_max+1)^2);
            N(idM) = N_subSet;
            N = tril(N) + tril(N)' + diag(diag(N));
            
            cLS = zeros(l_max+1);
            sLS = zeros(l_max+1);
            x = N\TN_subSet;
            i = 0;
            for az = 1 : l_max +1
                for m = 1 : az
                    i = i + 1;
                    cLS(az,m) = x(i);
                    if m > 1
                        i = i + 1;
                        sLS(az,m) = x(i);
                    end
                end
            end
            
            [ topo_map2 ] = sintesiGrid (phiGrid./180*pi, lambdaGrid./180*pi, cLS, sLS, 0, l_max, 0, l_max, 1, 1, 0, 0);
        end
        
        function plm = fplm(l, m, theta)
            % Computing Legendre polynomial
            %
            % INPUT
            %   l       degree
            %   m       order
            %   theta   theta angle [rad]
            %
            % SYNTAX
            %   plm = Core_Utils.fplm(l, m, decl)
            % 
            lMin = l(1);
            lMax = l(end);
            mMin = m(1);
            mMax = m(end);
            
            if (size(l,1)==1)
                l=l';
            end
            if (size(m,1)==1)
                m=m';
            end
            if (size(theta,1)==1)
                theta=theta';
            end
            
            %%
            r1 = zeros(mMax,1);
            
            % computing root 1 (bl)
            r1(1) = sqrt(3);
            i = (2:mMax);
            l = (1:lMax);
            
            r1(2:end) = sqrt((2*i+1)./(2*i));
            %%
            % Init P (result matrix)
            
            plm = zeros(mMax,length(m),length(theta));
            
            % Computing Pmm
            
            % Skip the calculous of the first l(1)-1 lines
            % Go from 0 to lMin-1
            Ptmp0 = ones(length(theta),1);
            for mfix = 1:mMin-1
                Ptmp1 = r1(mfix) * sin(theta).*Ptmp0;
                Ptmp0 = Ptmp1;
            end
            
            % for each m
            r2 = zeros(length(l)-1,1);
            r3 = zeros(length(l)-1,1);
            
            for mfix = mMin:mMax
                % Computing Pmm --------------------------------------------------
                Ptmp1 = r1(mfix) * sin(theta).*Ptmp0;
                Ptmp0 = Ptmp1;
                
                % Save in the results matrix the Pmm element
                plm(mfix, mfix-mMin+1, :) = Ptmp1(:);
                
                % Computing Plm --------------------------------------------------
                
                % get the row
                r = mfix+1;
                
                
                % computing root 2 (clm)
                r2 = sqrt(((2*l(r:end)+1).*(2*l(r:end)-1))./((l(r:end) + mfix).*(l(r:end) - mfix)));
                
                % computing root 3 (dlm)
                r3 = sqrt(((2*l(r:end)+1).*(l(r:end)+mfix-1).*(l(r:end)-mfix-1))./((2*l(r:end)-3).*(l(r:end)+mfix).*(l(r:end)-mfix)));
                
                Pl1 = Ptmp1;
                Pl2 = zeros(size(Ptmp1));
                for lfix = r:lMax
                    tmp = r2(lfix-r+1) * cos(theta).*Pl1 - r3(lfix-r+1).*Pl2;
                    plm(lfix,mfix-mMin+1, :) = tmp(:);
                    Pl2 = Pl1;
                    Pl1 = tmp(:);
                end
            end
            
            plm = plm(lMin:lMax,:,:);
        end

        %--------------------------------------------------------------------------
        %% TRIGONOMETRIC manipulators
        %--------------------------------------------------------------------------

        function r_angle = deg2rad(d_angle)
            % Convert degrees to radians
            %
            % SYNTAX
            %   r_angle =Core_Utils.deg2rad(d_angle)
            r_angle = d_angle .* (pi/180);
        end
        
        function d_angle = rad2deg(r_angle)
            % Convert degrees to radians
            %
            % SYNTAX
            %   r_angle =Core_Utils.deg2rad(d_angle)
            
            d_angle = r_angle .* (180/pi);
        end
                
        function [xe, ye] = ellipse(x, y, major, minor, az)
            % Get the coordinates of an ellipse
            %
            % SYNTAX
            %   [xe, ye] = Core_Utils.ellipse(x, y, major, minor, theta)
            center = [x, y];
            theta = 90 - az;

            % Generate x and y coordinates for the ellipse
            t = linspace(0, 2*pi);
            xe = center(1) + major*cos(t)*cosd(theta) - minor*sin(t)*sind(theta);
            ye = center(2) + major*cos(t)*sind(theta) + minor*sin(t)*cosd(theta);
        end

        function ph = plotEllipse(x, y, major, minor, az, varargin)
            % Plot the coordinates of an ellipse
            %
            % SYNTAX
            %   [xe, ye] = Core_Utils.plotEllipse(x, y, major, minor, theta)
            [xe, ye] = Core_Utils.ellipse(x, y, major, minor, az);
            if nargin < 6
                varargin = {};
            end
            if ~iscell(varargin)
                varargin = {varargin};
            end
            ph = plot(xe, ye, varargin{:});
        end

        function ph = patchEllipse(x, y, major, minor, az, color, varargin)
            % Patch the coordinates of an ellipse
            %
            % SYNTAX
            %   [xe, ye] = Core_Utils.plotEllipse(x, y, major, minor, theta)
            [xe, ye] = Core_Utils.ellipse(x, y, major, minor, az);
            if nargin < 7
                varargin = {};
            end
            if ~iscell(varargin)
                varargin = {varargin};
            end
            ph = patch(xe, ye, color, varargin{:});
        end

        function ph = m_plotEllipse(x, y, major, minor, az, varargin)
            % Plot the coordinates of an ellipse
            %
            % SYNTAX
            %   [xe, ye] = Core_Utils.plotEllipse(x, y, major, minor, theta)
            [xe, ye] = Core_Utils.ellipse(x, y, major, minor, az);
            if nargin < 6
                varargin = {};
            end
            if ~iscell(varargin)
                varargin = {varargin};
            end
            [xe, ye] = m_ll2xy(xe, ye);
            ph = plot(xe, ye, varargin{:});
        end

        function ph = m_patchEllipse(x, y, major, minor, az, color, varargin)
            % Patch the coordinates of an ellipse
            %
            % SYNTAX
            %   [xe, ye] = Core_Utils.plotEllipse(x, y, major, minor, theta)
            [xe, ye] = Core_Utils.ellipse(x, y, major, minor, az);
            if nargin < 7
                varargin = {};
            end
            if ~iscell(varargin)
                varargin = {varargin};
            end
            [xe, ye] = m_ll2xy(xe, ye);
            ph = patch(xe, ye, color, varargin{:});
        end
        
        %--------------------------------------------------------------------------
        %% ELEVATION MAPPING 
        %--------------------------------------------------------------------------

        function [elevation, flag] = getElevation(dlat_list, dlon_list, force_clean, sources)
            % Get elevation from an online DTM service:
            %
            % INPUT
            %   dlat_list   array of latitudes (degree)
            %   dlon_list   array of longitude (degree)
            %   force_clean flag T/F to clean the cache
            %   source      possible source: etopo1 (default) 'srtm30m', 'aster30m'
            %
            % SYNTAX
            %   [elevation, flag] = Core_Utils.getElevation(dlat_list, dlon_list, force_clean, sources)
            %
            % NOTE
            %   The function requires an internet connection

            persistent cache;

            if nargin < 4 || isempty(sources)
                sources = {'etopo1'};
            end
            % string to cell of strings
            if ischar(sources)
                sources = {sources};
            end

            if nargin < 3 || isempty(force_clean)
                force_clean = false;
            end


            if numel(dlat_list) > 100
                % I can only make 100 req at a time
                elevation = zeros(numel(dlat_list),1, 'single');
                flag = char(zeros(numel(dlat_list),1, 'uint8'));
                n_parts = ceil(numel(dlat_list)/100);

                for p = 1: n_parts
                    lim = [1 0] + min(100.*[p-1 p], numel(dlat_list));
                    
                    [elevation(lim(1):lim(2)), flag(lim(1):lim(2))] = Core_Utils.getElevation(dlat_list(lim(1):lim(2)), dlon_list(lim(1):lim(2)), force_clean, sources);
                    pause(1); % I need 1 seconds of pause between requests as a limit of the API
                end
            else
                % There are less than 100 coordinates to get
                elevation = [];
                flag = [];

                % Search in cache
                id_cache = zeros(numel(dlat_list),1, 'uint32');
                elevation = zeros(numel(dlat_list),2, 'single');
                for i = 1:numel(dlat_list)
                    lat = dlat_list(i);
                    lon = dlon_list(i);
                    if not(isempty(cache))
                        tmp = find((round(lat,3) == round(cache(:,1),4) & round(lon,4) == cache(:,2)), 1, 'first');
                        id_cache(i) = iif(isempty(tmp),0, tmp);
                    else
                        try
                            % Load cache from file
                            elevation_path = Core.getFilePath('elevation');
                            load(elevation_path, 'elevation_cache');
                            cache = elevation_cache;
                            tmp = find((round(lat,3) == round(cache(:,1),4) & round(lon,4) == cache(:,2)), 1, 'first');
                            id_cache(i) = iif(isempty(tmp),0, tmp);
                        catch ex
                            % Elevation_cache cache not found
                            cache = [];
                        end
                    end
                end

                if force_clean && any(id_cache > 0)
                    % clean the cache for the found points
                    id_tmp = id_cache(id_cache > 0);
                    id_ok = char(cache(id_tmp,4)) == sources{1}(2);
                    cache(id_tmp(~id_ok,:),:) = [];
                    id_cache(id_tmp(~id_ok)) = 0;
                end

                % Get what's already in cache
                if any(id_cache > 0)
                    elevation(id_cache > 0,:) = cache(id_cache(id_cache > 0), [3, 4]);
                end
                id_req = find(id_cache == 0); % Elevation to request

                try
                    if ~isempty(id_req)
                        status_ok = false;
                        cmd = '';
                        for i = 1:numel(sources)
                            if any(id_req)
                                source = sources{i};

                                loc_list = '';
                                for l = id_req(:)'
                                    loc_list = sprintf('%s|%.4f,%.4f', loc_list, dlat_list(l), dlon_list(l));
                                end
                                loc_list(1) = '';
                                cmd = sprintf('https://api.opentopodata.org/v1/%s?locations=%s', source, loc_list);
                                req = webread(cmd);
                                if req.status(1) == 'O'
                                    for r = 1:numel(req.results)
                                        elevation(id_req(r),:) = [single(req.results(r).elevation), single(source(2))];
                                    end
                                    if any(~isnan(elevation(id_req)))
                                        cache = [cache; round(dlat_list(id_req(~isnan(elevation(id_req)))),4), ...
                                            round(dlon_list(id_req(~isnan(elevation(id_req)))),4), ...
                                            elevation(id_req(~isnan(elevation(id_req))),:)];
                                        elevation_cache = cache; 
                                        elevation_path = Core.getFilePath('elevation');
                                        save(elevation_path, 'elevation_cache');                            
                                    end
                                end
                                id_req = id_req(isnan(elevation(id_req))); % if srtm is selected and results is nan => use aster
                            end
                        end
                    end
                catch ex
                    Core_Utils.printEx(ex);
                end
                % The 4th column of elevation is flag
                % flag is the second character of a requested source
                if ~isempty(elevation)
                    flag = char(elevation(:,2));
                    elevation = single(elevation(:,1));
                end
            end
        end

        %--------------------------------------------------------------------------
        %% COUNTRY MAPPING 
        %--------------------------------------------------------------------------
        
        function [country_iso3, country_iso2, country_name, iso_id] = getCountryISO3166(alpha_code, lon)
            % Get ISO codes of countries
            %
            % INPUT (variable)
            %   alpha_code2         2char code in standard ISO3166-1 alpha2
            %   alpha_code3         3char code in standard ISO3166-1 alpha3
            %   lat,lon             coordinates in degrees
            %
            % SYNTAX
            %   [country_iso3, country_iso2, country_name, iso_id] = Core_Utils.getCountryISO3166(alpha_code2);
            %   [country_iso3, country_iso2, country_name, iso_id] = Core_Utils.getCountryISO3166(alpha_code3);
            %   [country_iso3, country_iso2, country_name, iso_id] = Core_Utils.getCountryISO3166(lat, lon);
            
            persistent cache; % array with lat, lon, country_code
            
            if nargin == 1
                % getCountryISO3166(alpha_code)
                
                if isnumeric(alpha_code)
                    iso_id = alpha_code;
                elseif size(alpha_code,2) == 2
                    % alpha2 code
                    iso2code = Core_Utils.code2Char2Num(reshape([Core_Utils.ISO3166{:,1}], 2, size(Core_Utils.ISO3166, 1))');
                    for l = 1:size(alpha_code, 1)
                        iso_id(l) =  find(iso2code == Core_Utils.code2Char2Num(alpha_code(l,:)), 1, 'last');
                    end
                elseif size(alpha_code,2) == 3
                    % alpha3 code
                    iso2code = Core_Utils.code3Char2Num(reshape([Core_Utils.ISO3166{:,2}], 3, size(Core_Utils.ISO3166, 1))');
                    for l = 1:size(alpha_code, 1)
                        iso_id(l) =  find(iso2code == Core_Utils.code3Char2Num(alpha_code));
                    end
                end
                
                country_iso2 = reshape([Core_Utils.ISO3166{iso_id, 1}], 2, size(alpha_code, 1))';
                country_iso3 = reshape([Core_Utils.ISO3166{iso_id, 2}], 3, size(alpha_code, 1))';
                country_name = Core_Utils.ISO3166(iso_id, 3);
            elseif (nargin == 2)
                % getCountryISO3166(lat, lon)
                lat = alpha_code;
                
                country_iso3 = char(32*ones(numel(lat), 3, 'uint8'));
                country_iso2 = char(32*ones(numel(lat), 2, 'uint8'));
                country_name = cell(numel(lat), 1);
                for p = 1:numel(lat)
                    % search in cache
                    cached = false;
                    if not(isempty(cache))
                        id_cache = unique(cache((round(lat(p),3) == round(cache(:,1),3) & round(lon(p),3) == round(cache(:,2),3)), 3));
                        cached = any(id_cache);
                    else
                        try
                            % Load cache from file
                            nominatim_path = Core.getFilePath('nominatim');
                            load(nominatim_path, 'nominatim_cache');
                            cache = round(nominatim_cache,3);
                            id_cache = unique(cache((round(lat(p),3) == round(cache(:,1),3) & round(lon(p),3) == round(cache(:,2),3)), 3));
                            cached = any(id_cache);
                        catch ex
                            % nominatim cache not found
                            cache = [];
                        end
                    end
                    if cached
                        [country_iso3(p,:), country_iso2(p,:), country_name(p)] = Core_Utils.getCountryISO3166(id_cache);
                    else
                        if isnan(lat)
                            country_iso3(p,:) = Core_Utils.ISO3166{end,2};
                            country_iso2(p,:) = Core_Utils.ISO3166{end,1};
                            country_name(p,:) = Core_Utils.ISO3166(end,3);
                        else
                            [str_tmp] = Core_Utils.getLocationInfo(lat(p), lon(p));
                            
                            if isempty(str_tmp)
                                country_iso3(p,:) = Core_Utils.ISO3166{end,2};
                                country_iso2(p,:) = Core_Utils.ISO3166{end,1};
                                country_name(p,:) = Core_Utils.ISO3166(end,3);
                                str_tmp = struct('error', 'Unable to geocode');
                            end
                            if not(isempty(str_tmp))
                                flag_update_cache = false;
                                if isfield(str_tmp, 'error')
                                    Core.getLogger.addError(sprintf('%s (%g, %g)', str_tmp.error, lat(p), lon(p)));
                                    id_iso = size(Core_Utils.ISO3166,1); % world
                                    country_code = 'ERR';
                                    country_iso3(p,:) = 'ERR';
                                    country_iso2(p,:) = Core_Utils.ISO3166{id_iso,1};
                                    country_name(p) = Core_Utils.ISO3166(id_iso,3);
                                    flag_update_cache = true;
                                elseif isfield(str_tmp.address, 'country_code')
                                    country_code = upper(str_tmp.address.country_code);
                                    [country_iso3(p,:), country_iso2(p,:), country_name(p), id_iso] = Core_Utils.getCountryISO3166(country_code);
                                    flag_update_cache = true;
                                elseif isfield(str_tmp.address, 'man_made')
                                    Core.getLogger.addError(sprintf('Man made place "%s" (%g, %g)', str_tmp.address.man_made, lat(p), lon(p)));
                                    [country_iso3(p,:), country_iso2(p,:), country_name(p), id_iso] = Core_Utils.getCountryISO3166('OC');
                                    flag_update_cache = true;
                                else
                                    try
                                        tmp_addr = str_tmp.address.display_name;
                                    catch
                                        tmp_addr = 'unknown';
                                    end
                                    Core.getLogger.addError(sprintf('Man made place "%s" (%g, %g)', tmp_addr, lat(p), lon(p)));
                                    [country_iso3(p,:), country_iso2(p,:), country_name(p), id_iso] = Core_Utils.getCountryISO3166('WD');
                                    flag_update_cache = true;                                    
                                end
                                if flag_update_cache
                                    cache = round([cache; round(lat(p),3), round(lon(p),3), id_iso],3);
                                    fprintf('Adding to nominatim cache: %f, %f - %s\n', round(lat(p),3), round(lon(p),3), country_iso3(p,:))
                                    % load old nominatim, and remove repetitions
                                    % suboptimal intersect but faster
                                    try
                                        nominatim_path = Core.getFilePath('nominatim');
                                        load(nominatim_path, 'nominatim_cache');
                                    catch
                                        nominatim_cache = [];
                                    end
                                    nominatim_cache = round([nominatim_cache; cache],3); 
                                    [~, id_ok] = unique(nominatim_cache(:,1) + nominatim_cache(:,2)*1e8);
                                    nominatim_cache = nominatim_cache(id_ok,:);
                                    nominatim_path = Core.getFilePath('nominatim');
                                    save(nominatim_path, 'nominatim_cache');
                                    cache = nominatim_cache;
                                end
                            end
                       end
                    end
                end
            else
                error('Not the right amount of inputs');
            end
        end
        
        function [marker_name, id_num] = getMarkerV3(marker_name, lat, lon)
            % Return a markername in V3 format
            %
            % SYNTAX
            %   [marker_name, id_num] = Core_Utils.getMarkerV3(marker_name9ch)
            %   [marker_name, id_num] = Core_Utils.getMarkerV3(marker_name, coo)
            %   [marker_name, id_num] = Core_Utils.getMarkerV3(marker_name, lat, lon)
            
            if isnumeric(marker_name)
                % I suppose is a numeric marker_v3
                marker_name = Core_Utils.markerCode2MarkerName(marker_name);
            end
            
            marker_name = marker_name(1:min(9, numel(marker_name)));
            n_char = size(marker_name, 2);
            flag_v3 = Core_Utils.isMarkerV3(marker_name);
                
            if not(flag_v3)
                marker_name = strrep([marker_name repmat(' ',1,max(0,4 -length(marker_name)))], ' ', '_');%pad(marker_name, 4);
                marker_name = marker_name(1:4);
                if isa(lat, 'Coordinates')
                    [lat, lon] = lat.getMedianPos.getGeodetic();
                    lat = lat * 180/pi;
                    lon = lon * 180/pi;
                end
                country_iso3 = Core_Utils.getCountryISO3166(lat, lon);
                marker_name = [marker_name '00' country_iso3];
            end
            marker_name = [marker_name repmat(' ',1,max(0,9 -length(marker_name)))];%pad(marker_name, 9);
            if nargout == 2
                id_num = uint64(Core_Utils.code4Char2Num(marker_name(1:4)))*2^32 + uint64(str2double(marker_name(5:6))) * 2^24 + uint64(Core_Utils.code3Char2Num(marker_name(7:9)));
            end
        end
        
        function flag_v3 = isMarkerV3(marker_name)
            % return true if the marker is a valid RINEX3 name format
            %
            % SYNTAX
            %   flag_v3 = Core_Utils.isMarkerV3(marker_name)
            n_char = size(marker_name, 2);
            flag_v3 = true;
                
            if n_char == 9
                % Check if the format is V3
                % 5 and 6 char must be number between 0 and 9
                flag_v3 = flag_v3 && (marker_name(5) >= '0') && (marker_name(5) <= '9');
                flag_v3 = flag_v3 && (marker_name(6) >= '0') && (marker_name(6) <= '9');
                % Check if the current country code is valid
                iso2code = Core_Utils.code3Char2Num(reshape([Core_Utils.ISO3166{:,2}], 3, size(Core_Utils.ISO3166, 1))');
                flag_v3 = flag_v3 && not(isempty(find(iso2code == Core_Utils.code3Char2Num(marker_name(7:9)))));
            else
                flag_v3 = false;
            end
        end
        
        function marker_name = markerCode2MarkerName(id_num)
            % From a numeric marker get back the string format
            %
            % SYNTAX 
            %   marker_name = Core_Utils.markerCode2MarkerName(id_num)
            
            name = Core_Utils.num2Code4Char(bitand(id_num, uint64(sum(2.^(63:-1:32))))/2^32);
            ant_rec = reshape(sprintf('%02d', bitand(id_num, uint64(sum(2.^(24 + 0:7))))/2^24), numel(id_num), 2);
            country = Core_Utils.num2Code3Char(bitand(id_num, uint64(sum(2.^(23:-1:0)))));
            marker_name = [name, ant_rec, country];
        end
        
        function id_num = markerName2MarkerCode(marker_name)
            % From a numeric marker get back the string format
            %
            % SYNTAX 
            %   marker_name = Core_Utils.markerCode2MarkerName(id_num)
            
            [~, id_num] = Core_Utils.getMarkerV3(marker_name);
        end
        
        function info = searchPlace(country, city, cap)
            % Search for info (coordinates) o a place
            %
            % SINTAX
            %     info = Core_Utils.searchPlace(country, city, cap)   
            if nargin < 2
                city = '';
            end
            if nargin < 3
                cap = '';
            end
            request_url = ['https://nominatim.openstreetmap.org/search?' ...
                        iif(isempty(country), '', ['country=' country '&']), ...
                        iif(isempty(city), '', ['city=' city '&']), ...
                        iif(isempty(cap), '', ['postalcode=' cap '&']), ...
                        'format=json'];
                    
            info = webread(request_url,'Timeout',5);
        end
        
        %--------------------------------------------------------------------------
        %% MAP FUNCTIONS
        %--------------------------------------------------------------------------
        
        function [dlat, dlon, dlat_ext, dlon_ext, nwse] = getMapBoundaries(lat, lon, margin)
            % Given a set of points return the optimal boundaries for map visualization
            % INPUT
            %   dlat    latitude of the points [deg] [1 x n]
            %   dlon    latitude of the points [deg] [1 x n]
            %   margin  [lat lon] margin [deg]
            %
            % OUTPUT
            %   dlat_lim       latitude limits [deg] for the map viualization
            %   dlon_lim       longitude limits [deg] for the map viualization
            %   dlat_lim_ext   longitude limits [deg] for the map download and projections
            %   dlon_lim_ext   longitude limits [deg] for the map download and projections
            %
            % SYNTAX
            %  [dlat_lim, dlon_lim, dlat_lim_ext, dlon_lim_ext] = getMapBoundaries(dlat, dlon, <margin = 0>)
            if nargin < 3 || isempty(margin)
                margin = 0;
            end
            lat = noNaN(lat);
            lon = noNaN(lon);
            
            if numel(lon) == 1
                % If there is only one station
                dlon = minMax(lon) + [-0.05 0.05] + [-1 1] * margin(2);
                dlat = minMax(lat) + [-0.03 0.03] + [-1 1] * margin(1);
            else
                % If there are more stations add a border of 1/10 of the area among the points
                dlon = minMax(lon); dlon = dlon + [-1 1] * max(0.000015, margin(2) + (diff(dlon) / 2.5));
                dlat = minMax(lat); dlat = dlat + [-1 1] * max(0.00001, margin(1) + (diff(dlat) / 2.5));
            end
            if diff(dlon) <= 0
                dlon = dlon + [-1 1] * 0.00001;
            end
            if diff(dlat) <= 0
                dlat = dlat + [-1 1] * 0.00001;
            end
            dlat(1) = max(-90,dlat(1));
            dlat(2) = max(-90,dlat(2));
            if diff(dlon) > 360
                dlon = [-180.5 180.5];
            end
            % Set the image proportions to be max 3/4 (I don't like images that are too long)
            prop = 3/4;
            if diff(dlon) > diff(dlat) / cosd(mean(dlat))
                margin = max(diff(dlon) * prop * cosd(mean(dlat)), diff(dlat))/2;
                dlat = [dlat(1) + diff(dlat)/2 - margin, dlat(2) - diff(dlat)/2 + margin];
            else
                margin = max(diff(dlat) * prop / cosd(mean(dlat)), diff(dlon))/2;
                dlon = [dlon(1) + diff(dlon)/2 - margin, dlon(2) - diff(dlon)/2 + margin];
            end
            dlat = min(90, max(-90, dlat));
            
            nwse = [dlat(2), dlon(1), dlat(1), dlon(2)];
            dlon_ext = nwse([2 4]) + [-0.0001 0.0001];
            dlat_ext = max(-90, min(90, nwse([3 1]) + [-0.0001 0.0001]));
            
            if diff(dlat) < 0.0001  % if you make this too small drawing ticks in m_grid can screw up
                dlat = dlat + [-0.00005 +0.00005];
            end
            if diff(dlon) < 0.0001  % if you make this too small drawing ticks in m_grid can screw up
                dlon = dlon + [-0.00005 +0.00005];
            end
        end
        
        %--------------------------------------------------------------------------
        %% DATA manipulators
        %--------------------------------------------------------------------------
        function [data_cleaned, dt] = remCommonClock(data)
            dt = cumsum(nan2zero(robAdj(Core_Utils.diffAndPred(data))));
            data_cleaned = bsxfun(@minus, data, dt);
        end

        function [data_cleaned, id_ok, noise] = reduceCommonEffect(data, time, spline_base)
            % Reduce the common effect acting on all the columns of a dataset
            % e.g. tropospheric parameters from many stations in the same region
            %      contains a common effect due to satellites clock errors in PPP
            %
            % OUTPUT
            %   data_cleaned    data - common effect
            %   id_ok           column of data used to compute the common effect
            %   data_reduction  noise that have been reduced
            %
            % SYNTAX
            %   [data_cleaned, id_ok, data_reduction] = Core_Utils.reduceCommonEffect(data, time, spline_base)

            if nargin < 2
                spline_base = 3600;
            end
            data_cleaned = data;
            flag_debug = false;
            
            % Reduce each column by spline (hourly is the default), to keep the low degrees effects
            data_cleaned = data;
            id_ok = true(size(data,2),1);
            for r = 1:size(data,2)
                if any(data(:,r))
                    data_cleaned(:,r) = data(:,r) - splinerMat(time.getMatlabTime*86400, data(:,r), spline_base, 1e-9);
                else
                    id_ok(r) = false;
                end
            end

            % Compute STD for each column to find anomalous data
            std_data = std(data_cleaned, 'omitnan'); % Check for any anomalous data (too smooth => badly estimated)
            thr = mean(std_data(id_ok), 2, 'omitnan');
            id_ok = id_ok(:) & std_data(:) > thr/2 & std_data(:) < 3*thr ; % do not trust super low or high std

            % If I have at least 3 usable columns
            if sum(id_ok) > 3
                % Common factor is computed as a robust adjastment (Huber)
                noise = robAdj(data_cleaned(:,id_ok));
                % It should not be necessary, but it might be useful, reduce again the noise with the same spline used to reduce the data
                noise_splined = splinerMat(time.getMatlabTime*86400, noise, spline_base, 1e-9);

                if flag_debug
                    figure(100); clf;
                    plot(time.getMatlabTime, median(data_cleaned,2,'omitnan')); hold on;
                    plot(time.getMatlabTime, noise);
                    plot(time.getMatlabTime, noise_splined);
                end

                noise = noise - noise_splined;
                
                if flag_debug
                    figure(101); clf; subplot(2,1,1); plot(time.getMatlabTime, data); ax(1) = gca; subplot(2,1,2); plot(time.getMatlabTime, data - noise); ax(2) = gca; linkaxes(ax, 'x');
                    dockAllFigures;
                end
                
                data_cleaned = bsxfun(@minus, data, noise); %  compute "filtered" output
            end
        end

        function [t, data_set] = insertNan4Plots(t, data_set)
            % Insert a Nan in a regularly sampled dataset to make
            % plots interrupt continuous lines
            %
            % INPUT
            %   t      epoch of the data [matrix of column arrays]
            %   data   epoch of the data [matrix of column arrays]
            %
            % SYNTAX
            %   [t, data] = Core_Utils.insertNan4Plots(t, data)
              
            if isa(t, 'datetime')
                t = datenum(t);
                flag_datetime = true;
            else
                flag_datetime = false;
            end
            t = t(:);
            if size(t, 1) ~= size(data_set, 1)
                % data should be a column array
                data_set = data_set';
            end
            n_set = size(data_set, 2);
            dt = diff(t);
            rate = median(dt);
            id_in = find(dt > 1.5 * rate);
            id_in2 = id_in(:) + (1 : numel(id_in))';
            t2 = nan(numel(t)+numel(id_in2), 1);
            try
                id_ok = 1:numel(t2); id_ok(id_in2) = []; % find id in the new set of the old data
            catch
                keyboard
            end
            t2(id_ok) = t;
            data_set2 = nan(numel(t2), size(data_set,2));
            data_set2(id_ok, :) = data_set;
            data_set = data_set2;
            t = t2;
            %for x = numel(id_in) : -1 : 1
            %    t = [t(1 : id_in(x)); (t(id_in(x)) + 1.5 * rate); t((id_in(x)+1) : end)];
            %    data_set = [data_set(1 : id_in(x), :); nan(1, n_set); data_set((id_in(x)+1) : end, :)];
            %end
            if flag_datetime
                t = datetime(datevec(t));
            end
        end
        
        function [t, data_set] = insertZeros4Plots(t, data_set)
            % Insert a zero in a regularly sampled dataset to make
            % plots interrupt continuous lines
            %
            % INPUT
            %   t      epoch of the data [matrix of column arrays]
            %   data   epoch of the data [matrix of column arrays]
            %
            % SYNTAX
            %   [t, data] = Core_Utils.insertZeros4Plots(t, data)
                        
            t = t(:);
            if size(t, 1) ~= size(data_set, 1)
                % data should be a column array
                data_set = data_set';
            end
            n_set = size(data_set, 2);
            dt = diff(t);
            rate = median(dt);
            id_in = find(dt > 1.5 * rate);
            for x = numel(id_in) : -1 : 1
                t = [t(1 : id_in(x)); (t(id_in(x)) + 1e-9 * rate); (t(id_in(x) + 1) - 1e-9 * rate); t((id_in(x)+1) : end)];
                data_set = [data_set(1 : id_in(x), :); zeros(1, n_set); zeros(1, n_set); data_set((id_in(x)+1) : end, :)];
            end
        end
        
        function lh = plotSep(t, data, varargin)
            % Special wrapper to regular plot
            % Works on regularly sampled data
            % When there is a gap of data, it insert a nan value
            % to avoid linear interpolation of the data
            %
            % INPUT
            %   t           column array of epochs
            %   data        columns of data (could be a matrix)
            %   varagin     add other useful parameters of the plot
            %
            % SYNTAX
            %   lh = Core_Utils.plotSep(t, data, varagin);
            %
            % SEE ALSO
            %   plot
            if isa(t, 'handle')
                ax = t;
                t = data;
                data = varargin{1};
                varargin = varargin(2:end);
            else
                ax = gca;
            end
            
            flag_zeros = false;
            if nargin >= 3
                for i = 1 : numel(varargin)
                    if iscell(varargin)
                        if ischar(varargin{i}) && strcmp(varargin{i}, 'zeros')
                            flag_zeros = true;
                            varargin(i) = [];
                        end
                    else
                        if ischar(varargin) && strcmp(varargin, 'zeros')
                            flag_zeros = true;
                            varargin = [];
                        end
                    end
                end
            end
            try
                if numel(t) ~= numel(data) && numel(t) ~= size(data, 1)
                    % this means that there is only one data parameter and no time
                    varargin = [{data} varargin];
                    data = t;
                    if size(data, 1) == 1
                        % I want the data to be columnwise
                        data = data';
                    end
                    t = (1 : size(data, 1))';
                else
                    if size(data, 1) == 1
                        % I want the data to be columnwise
                        data = data';
                    end
                end
            catch
                % probably data is undefined
                data = t;
                if size(data, 1) == 1
                    % I want the data to be columnwise
                    data = data';
                end
                t = (1 : size(data, 1))';
            end
            
            for c = 1 : size(data, 2)
                if numel(data) == numel(t)
                    if flag_zeros
                        [t_col, data_col] = Core_Utils.insertZeros4Plots(t(:,c), data(:,c));
                    else
                        [t_col, data_col] = Core_Utils.insertNan4Plots(t(:,c), data(:,c));
                    end
                else
                    if flag_zeros
                    [t_col, data_col] = Core_Utils.insertZeros4Plots(t, data(:,c));
                    else
                    [t_col, data_col] = Core_Utils.insertNan4Plots(t, data(:,c));
                    end
                end
                if isempty(varargin)
                    lh = plot(ax, t_col, data_col);
                else
                    if iscell(varargin)
                        lh = plot(ax, t_col, data_col, varargin{:});
                    else
                        lh = plot(ax, t_col, data_col, varargin);
                    end
                end
                hold on;
            end
        end
        
        function ph_list = plotConfBand(t, data, err, color, alpha, varargin)
            % Special wrapper to regular plot
            % Works on regularly sampled data
            % When there is a gap of data, it closes the patch
            %
            % INPUT
            %   t           column array of epochs
            %   data        columns of data (could be a matrix)
            %   err         error level around data
            %   color       color of the patch
            %   alpha       alpha of the patch
            %   varagin     add other useful parameters of the plot
            %
            % SYNTAX
            %   lh = Core_Utils.plotConfBand(t, data, err, color, varargin);
            %
            % SEE ALSO

            if nargin < 3 || isempty(color)
                color = Core_UI.GRAY;
            end
            if nargin < 4 || isempty(alpha)
                alpha = 0.1;
            end
            ph_list = [];
            for c = 1:size(err,2)
                patch_lim = getFlagsLimits(~isnan(err(:,c)));
                for p = 1:size(patch_lim, 1)
                    x_patch = [t(patch_lim(p,1):patch_lim(p,2)); flipud(t(patch_lim(p,1):patch_lim(p,2)))]; % concatenate x values back-to-front
                    y_patch = [data(patch_lim(p,1):patch_lim(p,2),c) + err(patch_lim(p,1):patch_lim(p,2),c); flipud(data(patch_lim(p,1):patch_lim(p,2),c) - err(patch_lim(p,1):patch_lim(p,2),c))]; % concatenate y values for upper and lower bounds
                    ph_list = [ph_list patch(x_patch, y_patch, 'k', 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', alpha, 'HandleVisibility', 'off')];
                end
            end
        end

        function lh = patchSep(t, data, color, varargin)
            % Special wrapper to regular plot
            % Works on regularly sampled data
            % When there is a gap of data, it closes the patch
            %
            % INPUT
            %   t           column array of epochs
            %   data        columns of data (could be a matrix)
            %   color       color of the patch
            %   varagin     add other useful parameters of the plot
            %
            % SYNTAX
            %   lh = Core_Utils.patchSep(t, data, color, varargin);
            %
            % SEE ALSO
            %   plot
            
            
            if size(data, 1) == 1
                % I want the data to be columnwise
                data = data';
            end
            [t_col, data_col] = Core_Utils.insertZeros4Plots(t, data);
            data_col = nan2zero(data_col);
            if nargin > 3
                if iscell(varargin)
                    patch([t_col(:); t_col(end); t_col(1)], [data_col(:); 0; 0], color, varargin{:});
                else
                    patch([t_col(:); t_col(end); t_col(1)], [data_col(:); 0; 0], color, varargin);
                end
            else
                patch([t_col(:); t_col(end); t_col(1)], [data_col(:); 0; 0], color);
            end
        end
        
        function [h,a] = plotconf(x, y, sy, color)
            % plot data and confidence interval on top
            %
            % SYNTAX:
            %  [h,a] = Core_Utils.plotconf(x,y,sy,color)
            if not(isempty(y))
                h = plot(x,y,'color',color,'LineWidth',2);
                a = patch([x(:); flipud(x(:))],[y(:)-sy(:); flipud(y(:)+sy(:))],zeros(size([x(:); x(:)])),'FaceColor',color,'EdgeColor','none','FaceAlpha',0.2,'HandleVisibility','off');
            else
                h = [];
                y = 0;
                a = patch([x(:); flipud(x(:))],[y(:)-0*sy(:); flipud(y(:)+sy(:))],zeros(size([x(:); x(:)])),'FaceColor',color,'EdgeColor','none','FaceAlpha',0.2,'HandleVisibility','off');
            end
        end

        %--------------------------------------------------------------------------
        %% SHOW FUNCTIONS
        %--------------------------------------------------------------------------
        function fh_list = showVMFModel(vmf_file)
            % Show VMF file
            %
            % SYNTAX
            %  fh_list = Core_Utils.showVMFModel(vmf_file)
            
            try
                [txt, lim, has_cr] = Core_Utils.readTextFile(vmf_file);
                data_lim = str2num(txt((lim(6,1) + 22):lim(6,2))); % begin_lat begin_lon ah aw zhd zwd 
                
                table_dim = [diff(data_lim(1:2) / data_lim(5)) + 1, diff(data_lim(3:4) / data_lim(6)) + 1];
                if txt((lim(3,1))+25) == '1'
                    table_dim(2) = table_dim(2) - 1;
                end
                % Read the entire table in block
                data_table = reshape(str2num(strrep((serialize(txt(lim(8,1):lim(8 + table_dim(1) * table_dim(2) - 1, 2))))', newline, ' '))', 6, table_dim(1) * table_dim(2))';
                
                ah = zeros(table_dim(1), table_dim(2),'single');
                aw = zeros(table_dim(1), table_dim(2),'single');
                zhd = zeros(table_dim(1), table_dim(2),'single');
                zwd = zeros(table_dim(1), table_dim(2),'single');
                
                y = (data_table(:,1) - data_lim(1)) / data_lim(5) + 1;
                x = (data_table(:,2) - data_lim(3)) / data_lim(6) + 1;
                
                id = y + (x - 1) * table_dim(1);
                %clear x y;
                
                ah(id) = data_table(:,3);
                aw(id) = data_table(:,4);
                zhd(id) = data_table(:,5);
                zwd(id) = data_table(:,6);
                %clear id data_table txt lim;
                
                fh_list = figure;
                cmap = Cmap.gat2(128);
                subplot(2,2,1); imagesc(flipud(ah)*1e3); colorbar; colormap(cmap); title('ah'); caxis([1.1 1.3]);
                subplot(2,2,2); imagesc(flipud(aw)*1e3); colorbar; title('aw'); caxis([0.3 0.8]);
                subplot(2,2,3); imagesc(flipud(zhd)*1e2); colorbar; title('ZHD'); caxis([110 240]);
                subplot(2,2,4); imagesc(flipud(zwd)*1e2); colorbar; title('ZWD'); caxis([0 45]);

                fh_list = figure;
                cmap = Cmap.gat2(128);
                subplot(2,2,1); imagesc(flipud(zhd1)*1e2); colorbar; colormap(cmap); title('ZHD'); caxis([110 240]);
                subplot(2,2,2); imagesc(flipud(zwd1)*1e2); colorbar; title('ZWD'); caxis([0 45]);
                subplot(2,2,3); imagesc(flipud(zhd3)*1e2); colorbar; title('ZHD'); caxis([110 240]);
                subplot(2,2,4); imagesc(flipud(zwd3)*1e2); colorbar; title('ZWD'); caxis([0 45]);
                
                cmap = Cmap.gat2(128);
                figure; imagesc(flipud(zhd1 - zhd3)*1e2); colorbar; colormap(gat); title('ZHD diff'); %caxis([110 240]);
                figure; imagesc(flipud(zwd1 - zwd3)*1e2); colorbar; title('ZWD diff'); %caxis([0 45]);
                
            catch ex
                Core_Utils.printEx(ex);
                Core.getLogger.addWarning(sprintf('Malformed VMF file "%s"', vmf_file));
            end
        end

        function setFigureName(fh, fig_name)
            % Set the name of the figure
            %
            % SYNTAX
            %   setFigureName(fh, fig_name)
            
            fh.NumberTitle = 'off'; fh.Name = [num2str(fh.Number, '%03d:') sprintf(' %s', fig_name)];
        end

        %--------------------------------------------------------------------------
        %% MEMORY
        %--------------------------------------------------------------------------
        
        function total_mem = getMem(obj, indent_lev)
            % Get all properties
            props = properties(obj);
            if nargin == 1
                indent_lev = 0;
            end
            if isempty(props)
                cur_prop = obj;
                s = whos('cur_prop');
                total_mem = s.bytes;
            else
                total_mem = 0;
                % Loop properties
                for p = 1 : length(props)
                    % Make shallow copy
                    cur_prop = obj.(props{p});  %#ok<*NASGU>
                    % Get info struct for current property
                    if isobject(cur_prop)
                        cur_mem = Core_Utils.getMem(cur_prop, indent_lev + 1);
                        fprintf('%s  ^ %s %s\n', char(32 * ones(1, indent_lev * 1, 'uint8')), byte2human(cur_mem), props{p});
                    else
                        s = whos('cur_prop');
                        % Add to total memory consumption
                        cur_mem = s.bytes;
                        % fprintf('%s |- %g bytes %s\n', char(32 * ones(1, indent_lev * 1, 'uint8')), cur_mem, props{p});
                    end
                    total_mem = total_mem + cur_mem;
                end
            end
            if indent_lev == 0
                fprintf('TOTAL occupation: %s\n', byte2human(total_mem));
            end

            function str = byte2human(n_bytes)
                if n_bytes > 1e9
                    str = sprintf('%g GB', n_bytes / 1e9);
                elseif n_bytes > 1e6
                    str = sprintf('%g MB', n_bytes / 1e6);
                elseif n_bytes > 1e3
                    str = sprintf('%g KB', n_bytes / 1e3);
                else
                    str = sprintf('%g bytes', n_bytes);
                end
            end
        end
        
        %--------------------------------------------------------------------------
        %% FILE MANAGEMENT
        %--------------------------------------------------------------------------
        
        function [txt, lim, has_cr] = readTextFile(file_path, min_line_len)
            % Read a full file as txt
            % 
            % INPUT: 
            %   file_path      full path to the text file
            %   min_line_len   lines shorter than this value are not considered (default = 3)
            %
            % OUTPUT:
            %   txt            full file as a string
            %   lim            limits of the lines [character of start stop, line width]
            %   has_cr         has carriage return character
            %
            % SYNTAX:
            %  [txt, lim, has_cr] = Core_Utils.readTextFile(file_path, min_line_len)
            
            fid = fopen(file_path, 'rt');
            if fid <= 0
                Core.getLogger.addError(sprintf('"%s" cannot be read', file_path));
                txt = [];
                lim = [];
                has_cr = false;
            else
                if nargin == 1
                    min_line_len = 3;
                end
                txt = fread(fid,'*char')';
                % try to see if carriage return is present in the file (Windows stupid standard)
                % On Windows file lines ends with char(13) char(10)
                % instead of just using char(10)
                if ~isempty(find(txt(1:min(1000,numel(txt))) == 13, 1, 'first'))
                    has_cr = true;  % The file has carriage return - I hate you Bill!
                else
                    has_cr = false;  % The file is UNIX standard
                end
                % txt = txt(txt ~= 13);  % remove carriage return - I hate you Bill!
                fclose(fid);
                
                % get new line separators
                nl = regexp(txt, '\n')';
                if (isempty(nl))
                    lim = [];
                else
                    if nl(end) <  (numel(txt) - double(has_cr))
                        nl = [nl; numel(txt)];
                    end
                    lim = [[1; nl(1 : end - 1) + 1] (nl - 1 - double(has_cr))];
                    lim = [lim lim(:,2) - lim(:,1) + 1];
                    l = 1;
                    while l <= size(lim,1)
                        if lim(l,3) < min_line_len
                            lim(l,:) = [];
                        else
                            l = l + 1;
                        end
                    end
                    
                    % removing empty lines at end of file
                    if size(lim,1) > 1
                        while lim(end,3)  < 1
                            lim(end,:) = [];
                        end
                    end
                end
            end
        end
        
        %--------------------------------------------------------------------------
        %% DEBUGGING FUNCTIONS
        %--------------------------------------------------------------------------
        
        function printEx(ex)
            % Print exception to screen
            %
            % SYNTAX:
            %  Core_Utils.printEx(ex)
            
            msg = sprintf(['goGPS encountered a problem, please open an issue on GitHub posting\n' ...
                    ' the following lines together with a copy of the full log' ]);
            str = sprintf('\n---------------------------------------------------------------------\n %s', msg);
            str = sprintf('%s\n---------------------------------------------------------------------\n MESSAGE: %s\n---------------------------------------------------------------------\n\n', str, ex.message);
            
            for i=1:numel(ex.stack)
                str = sprintf('%s  file: "%s"\n  line: %d\n  fun: %s\n\n', str, ex.stack(i).file, ex.stack(i).line, ex.stack(i).name);
            end
            fprintf(str);
            Core.getLogger.addMonoMessage(str);
            % keyboard
        end
                
        %--------------------------------------------------------------------------
        %% UTILITIES FUNCTIONS
        %--------------------------------------------------------------------------

        function [p] = getPositiveEntries(A)
            % Return what entry are positive in a matrix
            % 
            % OUTPUT
            %   p   logical array when an entry is positive
            %
            % SYNTAX
            %   [p] = getPositiveEntries(A);
            p = zeros(size(A,1),1); 
            for i = 1:size(A,1)
                tmp = A(1:i, 1:i); 
                tmp(:, p > 0) = []; 
                tmp(p > 0, :) = []; 
                [~, p(i)] = chol(tmp); 
            end
            p = p == 0;
        end
        
        function uuid = getUUID()
            % Get a unique identifier (char)
            %
            % SYNTAX
            %   Core_Utils.getUUID();
            
            uuid = strrep(char(java.util.UUID.randomUUID),'-','');
        end
        
        function sftp_par = getParSFTP(sftp_name)
            % Get connection parameter for a SFTP server 
            % read them from credentials.txt
            %
            % SYNTAX
            %   sftp_par = Core_Utils.getParSFTP(sftp_name)
            sftp_par = struct('name', sftp_name, ...
                'remote_host', '', ...
                'type', '', ...
                'port', 22, ...
                'username', '', ...
                'password', '', ...
                'folder', '');

            % read credentials
            credentials_path = Core.getFilePath('credentials');
            if exist(credentials_path, 'file') == 0
                Core.getLogger.addError('The "credentials.txt" file is missing.\nIt will be created as empty from "credentials.example.txt" in Breva folder');
                try
                    credentials_default_path = [Core.getInstallDir filesep 'credentials.example.txt'];
                    copyfile(credentials_default_path, credentials_path);
                    credentials = Ini_Manager(credentials_path);
                    credentials.readFile();
                catch
                    credentials = Ini_Manager();
                end
            else
                credentials = Ini_Manager(credentials_path);
                credentials.readFile();
            end

            sftp_par.remote_host = credentials.getData(['rs-' sftp_name],'remote_host');
            sftp_par.port = credentials.getData(['rs-' sftp_name],'port');
            sftp_par.username = credentials.getData(['rs-' sftp_name],'username');
            sftp_par.password = credentials.getData(['rs-' sftp_name],'password');
            sftp_par.folder = credentials.getData(['rs-' sftp_name],'folder');
            sftp_par.type = credentials.getData(['rs-' sftp_name],'type');
        end

        function pushFile(rem_loc_id, file_path, rem_folder)
            % Push a file on a remote location
            % At the moment it only support sftp servers
            %
            % INPUT
            %   rem_loc_id    name of the entry in credentials.txt with the parameters of the remote server
            %   file_path     path to the file to push on sftp
            %   rem_folder    folder (it must be existing) on the remote server where to push the file,
            %                 default = folder parameter as read from the configuration in credetials.txt
            %
            % SYNTAX
            %   Core_Utils.pushFile(rem_loc_id, file_path, <rem_folder>)
            
            % check if there is a folder attached to the rem_loc_id           
            rem_folder_postfix = regexp(rem_loc_id, "/.*", 'match', 'once');
            if ~isempty(rem_folder_postfix)
                rem_loc_id(end-length(rem_folder_postfix)+1 : end) = [];
            end
            % get the parameters for the connection
            sftp_par = Core_Utils.getParSFTP(rem_loc_id);
            
            % set the remote path as relative to the parameter in the location ini file
            if nargin < 3 || isempty(rem_folder)
                rem_folder = sftp_par.folder;
            end
            if isempty(rem_folder)
                rem_folder = '/';
            end
            if ~isempty(rem_folder_postfix)
                rem_folder = fullfile(rem_folder, rem_folder_postfix);
            end
            if isempty(sftp_par.remote_host)
                Core.getLogger.addError(sprintf('Remote location "%s" not found', rem_loc_id));
            else
                if isstring(file_path)
                    file_path = {file_path};
                end
                last_send = 0;
                try
                    sh = sftp([sftp_par.remote_host ':' num2str(sftp_par.port)], sftp_par.username, "password", sftp_par.password);
                    % Try to connect
                    try
                        fnp = File_Name_Processor; 
                        core = Core.getCurrentCore; 
                        rem_folder = fnp.dateKeyRep(rem_folder, core.getSessionCentralTime, core.getCurSession);
                        folder_list = regexp(fullfile(rem_folder,'/'),'(?<=/)[^/]*(?=/)', 'match');
                        if ~isempty(folder_list)
                            path = '/';
                            for f = 1:numel(folder_list)
                                path = fullfile(path, folder_list{f});
                                try
                                    mkdir(sh, path);
                                catch ex
                                    % the path exist / no permission
                                end
                            end
                        end
                        cd(sh, rem_folder);
                        for f = 1:numel(file_path)
                            Core.getLogger.addMessage(sprintf('Pushing "%s" to "%s://%s%s"\n', file_path{f}, sftp_par.type, [sftp_par.remote_host ':' num2str(sftp_par.port)], rem_folder));
                            mput(sh, file_path{f});
                            last_send = last_send + 1;
                            Core.getLogger.addStatusOk(sprintf('"%s" sent!', file_path{f}));
                        end
                    catch ex
                        Core_Utils.printEx(ex);
                        for f = (last_send+1):numel(file_path)
                            Core.getLogger.addError(sprintf('Something went wrong "%s" not sent!', file_path{f}));
                        end
                        return
                    end
                catch
                    Core_Utils.printEx(ex);
                    for f = (last_send+1):numel(file_path)
                        Core.getLogger.addError(sprintf('Something went wrong "%s" not sent!\nCheck your server parameters', file_path{f}));
                    end
                    return
                end
                close(sh);
            end
        end

        function coords = selectPolygon(lat, lon)
            % selectPolygon - Selects a polygon from a map using mouse clicks.
            %
            % This function allows users to select a polygon on a map. Users can
            % add vertices to the polygon by right-clicking and can drag existing
            % vertices by left-clicking. A vertex changes color when it is near the cursor
            % and can be dragged. Instructions and keyboard shortcuts are provided for ease of use.
            %
            % INPUTS:
            % - lat, lon: Latitude and longitude to center the map (optional).
            %
            % OUTPUTS:
            % - coords: A matrix containing the coordinates of the selected polygon.

            % Check the number of arguments to set the center of the map
            if nargin == 1 
                [lat, lon] = lat.getMedianPos.getGeodetic();
                lat = lat * 180 / pi;
                lon = lon * 180 / pi;                
            end

            % Create figure and set map limits
            fh = figure;
            xlim([lon-0.002, lon+0.002]);
            ylim([lat-0.002, lat+0.002]);
            Core_UI.beautifyFig(fh);
            gm = addMap('alpha', 0.95, 'provider', 'satellite');
            plot(lon, lat, '^', 'color', Core_UI.ORANGE, 'markersize', 10, 'linewidth', 5); hold on;
            instructions = sprintf('Right-Click to add points.\nLeft-Click near point to drag');
            title(sprintf('Design a polygon\n\\fontsize{10} %s\\fontsize{5} \n',instructions));

            % Instructions
            
            % Initialize variables for polygon and UI
            coords = [];
            h_line = line(nan, nan, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '-.', 'Marker', '.', 'MarkerSize', 30, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'Parent', gca);
            h_patch = patch(nan, nan, 'b', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'FaceAlpha', 0.3, 'Parent', gca);
            h_dragged_point = plot(nan, nan, 'o', 'color', 'red', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'w', 'MarkerSize', 12, 'LineWidth', 2); % Dragged point
            uistack(h_dragged_point, 'top'); % Bring the dragged point to the front
            txt = uicontrol('Style', 'edit', 'Position', [60 5 400 20]);
            fig_width = fh.Position(3);
            txt.Position(3) = fig_width - txt.Position(1) - 10;
            drawnow;

            % Set up listeners for mouse interactions and window resize
            set(gm, 'ButtonDownFcn', @mouseClickFcn);
            set(h_patch, 'ButtonDownFcn', @mouseClickFcn);
            set(fh, 'WindowButtonUpFcn', @mouseReleaseFcn);
            set(fh, 'KeyPressFcn', @keyPressFcn);
            set(fh, 'SizeChangedFcn', @resizeUI);

            % Initialize dragging-related variables
            is_dragging = false;
            dragged_vertex_index = [];

            % Add UI controls for undo functionality
            undu = uicontrol('Style', 'pushbutton', 'String', 'Undo', 'Position', [20, 5, 40, 20], 'Callback', @undoCallback);

            % Mouse click function to handle adding and dragging points
            function mouseClickFcn(~, event)
                pt = get(gca, 'CurrentPoint');
                pt = [pt(1,1), pt(1,2)];

                [is_near, index] = isNearVertex(pt, coords);
                if event.Button == 1 && is_near % Left button for dragging
                    is_dragging = true;
                    dragged_vertex_index = index;
                    set(h_dragged_point, 'XData', coords(index,1), 'YData', coords(index,2), 'color', 'r');
                    set(gcf, 'WindowButtonMotionFcn', @mouseDragFcn);
                elseif event.Button == 3 % Right button to add new vertex
                    coords(end+1,:) = pt;
                    updateLine();
                end
            end

            % Function to check if a point is near any vertex of the polygon
            function [is_near, index] = isNearVertex(point, vertices)
                if isempty(vertices)
                    is_near = false;
                    index = 0;
                else
                    tolerance = 0.001; % Define a tolerance for proximity
                    distances = sqrt(sum((vertices - point).^2, 2));
                    [min_distance, index] = min(distances);
                    is_near = min_distance < tolerance;
                end
            end

            % Function to handle dragging of a vertex
            function mouseDragFcn(~, ~)
                if is_dragging
                    pt = get(gca, 'CurrentPoint');
                    pt = [pt(1,1), pt(1,2)];
                    coords(dragged_vertex_index, :) = pt;
                    set(h_dragged_point, 'XData', pt(1), 'YData', pt(2));
                    updateVisuals();
                end
            end

            % Function to handle release of mouse button after dragging
            function mouseReleaseFcn(~, ~)
                if is_dragging
                    is_dragging = false;
                    set(gcf, 'WindowButtonMotionFcn', '');
                    set(h_dragged_point, 'XData', nan, 'YData', nan, 'color', 'green'); % Reset the dragged point
                    updateLine();
                end
            end

            % Function to update the visuals of the line and patch
            function updateVisuals()
                                set(h_line, 'XData', coords([1:end,1],1), 'YData', coords([1:end,1],2));
                set(h_patch, 'XData', coords([1:end,1],1), 'YData', coords([1:end,1],2));
                end

            % Function to update the line, patch, and text field
            function updateLine()
                updateVisuals();
                set(txt, 'String', sprintf('poly_mask = %s', mat2str(coords, 9)));
            end

            % Undo callback to remove last point
            function undoCallback(~, ~)
                if ~isempty(coords)
                    coords(end,:) = [];
                    updateLine();
                end
            end

            % Function to handle keyboard shortcuts
            function keyPressFcn(~, event)
                if strcmp(event.Key, 'z') && any(strcmp(event.Modifier, 'control'))
                    undoCallback();
                end
            end

            % Function to resize UI elements
            function resizeUI(~, ~)
                % Adjust positions and sizes of UI elements based on the figure size
                new_width = fh.Position(3);
                txt.Position(3) = new_width - txt.Position(1) - 10;
            end

            % Wait for user to press enter or close the figure
            uiwait(fh);
        end
        
        function in_poly_lid = inPolygon(poly_xy, x, y)
            % This function returns a logical index for the points that fall within a polygon.
            % Using cartesian approximation
            %
            % INPUTS:
            % - coords: an n x 2 matrix containing the longitude and latitude coordinates of the polygon vertices
            % - x:      an m x 1 array containing the longitude coordinates of the points
            % - y:      an m x 1 array containing the latitude coordinates of the points
            %
            % OUTPUTS:
            % - in_poly_lid: an m x 1 logical array indicating whether each point falls within the polygon
            %
            % SYNTAX
            %   in_poly_lid = Core_Utils.inPolygon(poly_xy, x, y)
            
            in_poly_lid = inpolygon(x, y, poly_xy(:,1), poly_xy(:,2));
        end

        %--------------------------------------------------------------------------
        %% EXPORT FUNCTIONS
        %--------------------------------------------------------------------------
        
        function tropo = tropoAggregator(search_dir, flag_figure, flag_save)
            % Aggregate .mat tropo export from different receiver into a unique struct
            %
            % SYNTAX
            %   tropo = tropoAggregator(search_dir, flag_figure, flag_save);
            tropo = [];
            if nargin < 2 || isempty(flag_save)
                flag_save = true;
            end
            if nargin < 3 || isempty(flag_save)
                flag_save = false;
            else
                if ischar(flag_save)
                    export_path = flag_save;
                    flag_save = true;
                else
                    export_path = fullfile(Core.getState.getOutDir, [Core.getState.getPrjName, '_', GPS_Time.now.toString('yyyymmdd'), '.mat']);
                end                
            end
            min_num_epochs = 0;
            tropo = struct();
            files = dir(fullfile(search_dir, '*.mat'));
            lat_lim = [-90 90];
            lon_lim = [-180 180];
            if flag_figure
                figure;
            end
            for i = 1:length(files)
                fname = files(i).name;
                if strfind(fname,'.mat') & (length(fname) > 4)
                    Core.getLogger.addMessage(sprintf('Aggregating %s', fname));
                    tropo_data = load([ search_dir '/' fname]);
                    ztd_ok = false;% p(1) + p(2) * tropo_data.h_ortho + p(3) * tropo_data.h_ortho.^2;
                    d_ztd = 0;%diff(tropo_data.ztd)./diff(tropo_data.utc_time)*(30/86400);
                    mv_avg = ones(1,20)/20;
                    md_ztd = filter(mv_avg,1,d_ztd);
                    lon = tropo_data.lon;
                    lon(lon > 180) = lon(lon > 180)-360;

                    if false && sum(abs(tropo_data.ztd - ztd_ok) > 0.3) > 1 || sum(isnan(tropo_data.ztd)) > 1 || length(tropo_data.ztd) < min_num_epochs || max(abs(md_ztd)) > 0.005
                        % disp([' Outlier!! ' fname(1:4)])
                        if max(abs(md_ztd)) > 0.005 & sum(abs(tropo_data.ztd - ztd_ok) > 0.3) == 0 || sum(abs(tropo_data.ztd - ztd_ok) > 0.3) > 1
                            disp([' Outlier!! ' fname(1:4)])
                        else
                            disp([' Empty or partly empty!! ' fname(1:4)])
                            disp( [sum(isnan(tropo_data.ztd)) length(tropo_data.ztd)])
                        end
                    elseif tropo_data.lat > lat_lim(1) & tropo_data.lat < lat_lim(2) & lon > lon_lim(1) & lon < lon_lim(2)
                        if flag_figure
                            scatter(lon,tropo_data.lat);
                            hold on;
                        end
                        if isfield(tropo,fname(1:4))
                            if abs(tropo.(fname(1:4)).h_ellips - tropo_data.h_ellips) > 0.1
                                sprintf('Warning dh = %f',(tropo.(fname(1:4)).h_ellips - tropo_data.h_ellips) );
                            end
                            tropo.(fname(1:4)).ztd = [tropo.(fname(1:4)).ztd ; tropo_data.ztd];
                            tropo.(fname(1:4)).utc_time = [tropo.(fname(1:4)).utc_time ; tropo_data.utc_time];
                        else
                            tropo.(fname(1:4)) = tropo_data;
                        end
                    end
                end
            end
            if ~isempty(tropo) && ~isempty(files)
                Core.getLogger.addMarkedMessage(sprintf('Exporting to %s', export_path));
                save(export_path, 'tropo');
            end
            if isempty(files)
                Core.getLogger.addError(sprintf('No data .mat found in "%s"', search_dir));
            end
        end

        function exportTropo()
        end

        function exportCurFig(out_path, mode)
            % Eport Current Figure
            %
            % EXAMPLE
            %   Core_Utilis.exportCurFig(fullfile('/Users/Andrea/Library/Mobile Documents/com~apple~CloudDocs/MyDocuments/GIMS/title.png'))
            if nargin == 2
                Core_Utils.exportFig(gcf, out_path, mode);
            else
                Core_Utils.exportFig(gcf, out_path);
            end
        end
        
        function exportFig(fh, out_path, mode, flag_transparent)
            % Export Figure
            %
            % SYNTAX
            %   Core_Utils.exportFig(<fh>, out_path)
            %   Core_Utils.exportFig(out_path)
            %
            % EXAMPLE
            %   Core_Utilis.exportFig(fh, fullfile('/Users/Andrea/Library/Mobile Documents/com~apple~CloudDocs/MyDocuments/GIMS/title.png'))
            if nargin <= 1 
                out_path = fh;
                fh = gcf;
            end            
            %fh.WindowStyle = 'normal'; export_fig(fh, out_path, '-transparent', '-r150'); fh.WindowStyle = 'docked';
            ws_bk = fh.WindowStyle;
            fh.WindowStyle = 'normal';
            if nargin >= 3 && ~isempty(mode)
                Core_UI.beautifyFig(fh, mode);
            end
            if (nargin < 4) || isempty(flag_transparent)
                as = App_Settings.getInstance;
                as.reload;
                flag_transparent = as.isExportTransparent();
            end
            col = fh.Color;
            Logger.getInstance.addMessage(sprintf('Exporting to "%s"', out_path));
            box = findall(fh, 'type', 'uicontainer');
            if isempty(box)
                if flag_transparent
                    export_fig(fh, out_path, '-transparent', '-r150','-background', 'none');
                else
                    if fh.Visible
                        export_fig(fh, out_path, '-r150');
                    else
                        if ~isdeployed
                            original_position = get(fh, 'Position');
                            set(fh, 'Position', original_position + [+5e5 +5e5 0 0]); % Resize the figure back to its original size
                            fh.Visible = 'on';
                        end
                        export_fig(fh, out_path, '-r150');
                        if ~isdeployed
                            fh.Visible = 'off';
                            set(fh, 'Position', original_position);
                        end
                    end
                end
            else
                % Special tricks in case of figure containing boxes
                
                % Use saveas instead of export_fig
                [~, ~, ext] = fileparts(out_path);
                if strcmp(ext, '.png')
                    bg_color = [1 255 1]; % Use green screen
                else
                    bg_color = [255 255 255];
                end
                % fallback
                bg_box = {};
                for b = 1: numel(box)
                    try
                        bg_box{b} = box(b).BackgroundColor;
                        box(b).BackgroundColor = bg_color/255;
                    catch ex
                    end
                end
                saveas(fh, out_path);
                for b = 1: numel(box)
                    try                        
                        box(b).BackgroundColor = bg_box{b};
                    catch ex
                    end
                end
                
                % saveas does not have transparency management
                if strcmp(ext, '.png') && flag_transparent
                    % Tricks are just for PNG file type

                    % read the just saved image 
                    [im_out] = imread(out_path);
                    % convert it in hue saturation value
                    im_hsv = rgb2hsv(im_out);
                    bg_hsv = rgb2hsv(bg_color/255);
                    alpha_mask = (im_hsv(:,:,1) ~= bg_hsv(1));
                    color_mask = (~(im_out(:,:,1) == bg_color(1) & im_out(:,:,2) == bg_color(2) & im_out(:,:,3) == bg_color(3)));
                    saturation = im_hsv(:,:,2);
                    value = im_hsv(:,:,3);
                    hue = im_hsv(:,:,1);
                    saturation(alpha_mask ~= 1) = 0;
                    hue(alpha_mask ~= 1) = 1;
                    alpha_mask = double(alpha_mask);
                    
                    % get alpha value for the value component of the image
                    alpha_mask(logical(color_mask - alpha_mask)) = double(1-value(logical(color_mask - alpha_mask)));
                    
                    im_hsv(:,:,1) = hue;
                    im_hsv(:,:,2) = saturation;
                    im_hsv(:,:,3) = value;
                    im_out = hsv2rgb(im_hsv);
                    
                    imwrite(im_out, out_path, 'alpha', alpha_mask);
                end               
            end
                
            fh.WindowStyle = ws_bk;
            fh.Color = col;
        end
        
        function export2Plotly(fh, flag_offline, new_rate)
            % Export figure on plotly managing time axes
            %
            % SYNTAX
            %   export2Plotly(fig_handle, flag_offline, new_rate);
            %
            % NOTE 
            %   to use this function you must install plotly API:
            %       https://plot.ly/matlab/getting-started/#installation
            
            
            if nargin < 3
                new_rate = [];
            end
            
            if nargin < 2 || isempty(flag_offline)
                flag_offline = true;
            end
            
            if nargin > 1 && ~isempty(fh)
                figure(fh);
            else
                fh = gcf;
            end
            
            % Prepare times for plotly
            ax = findall(fh, 'type', 'axes');
            last = struct();
            last.mode = {};
            for i = 1 : numel(ax)
                last.mode{i} = ax(i).XTickLabelMode;
                ax(i).XTickLabelMode = 'auto';
            end
            
            line = findall(fh, 'type', 'line');
            last.flag_date = [];
            for i = 1 : numel(line)
                if ~isempty(new_rate)
                    time = line(i).XData * 86400;
                    rate = median(round(diff(time), 3));
                    [~, id_ok] = ismember((unique(round(time ./ new_rate)) .* new_rate), (round(time ./ (2*rate)) .* (2*rate)));
                    line(i).XData = line(i).XData(noZero(id_ok));
                    line(i).YData = line(i).YData(noZero(id_ok));
                end
                
                %if all(line(i).XData > datenum('1970-01-01') & line(i).XData < datenum('2070-01-01'))
                    %line(i).XData = convertDate(line(i).XData);
                    last.flag_date(i) = true;
                %else
                %    last.flag_date(i) = false;
                %end
            end            
            
            % Export to plotly
            try
                fig_name = fh.UserData.fig_name;
            catch
                fig_name = 'temp_unknown';
            end
            plotlyfig = fig2plotly(gcf, 'filename', fig_name, 'offline', flag_offline);

            if any(last.flag_date)
                for i = 1 : numel(ax)
                    plotlyfig.layout.(sprintf('xaxis%d', i)).type = 'date';
                end
                plotlyfig.PlotOptions.FileOpt = 'overwrite';
                plotly(plotlyfig);
            end
            
            for i = 1 : numel(line)
                if last.flag_date(i)
                    funToMatTime = @(date) (date/(1000*60*60*24) + datenum(1969,12,31,19,00,00));
                    line(i).XData = funToMatTime(line(i).XData);
                end
            end
            
            % Restore time for matlab
            for i = 1 : numel(ax)
                ax(i).XTickLabelMode = last.mode{i};
                setTimeTicks(ax(i));
            end
            
            Core_UI.beautifyFig(fh);
        end
                
        %--------------------------------------------------------------------------
        %% OTHER FUNCTIONS
        %--------------------------------------------------------------------------

        function ftp_server = ftpTimeout(ip, time_out_seconds, user_name, passwd)
            % Call FTP with timeout
            %
            % INPUT
            %   ip                  ip address
            %   time_out_seconds    time out [in seconds]
            %   user_name           (optional username)
            %   passwd              (optional password)
            %
            % SYNTAX
            %   ftp_server = ftpTimeout(ip, time_out_seconds, <user_name>, <passwd>)
            %

            import org.apache.commons.net.ftp.FTPClient;
            import java.net.SocketTimeoutException;

            time_out_seconds = time_out_seconds * 1000;
            % Attempt to connect to the FTP server using java
            try
                ftp_client = FTPClient();
                ftp_client.setConnectTimeout(time_out_seconds);
                ftp_client.connect(strrep(ip,':21', ''));
                ftp_client.logout();
                ftp_client.disconnect();
                
                if nargin == 4
                    ftp_server = ftp(ip, user_name, passwd);
                else
                    ftp_server = ftp(ip);
                end
                % Use ftp_server for your operations here
            catch e
                ftp_server = [];
            end
        end

        function str = exportNominatimCache(file_path)
            % Export in text file the nominatim cache
            % Load cache from file
            nominatim_path = Core.getFilePath('nominatim');
            load(nominatim_path, 'nominatim_cache');
            str = sprintf('   Lat    |    Lon    | CODE\n');
            for i = 1 : size(nominatim_cache,1)
                str = sprintf('%s %8.4f | %9.4f | %3s\n', str, nominatim_cache(i,1), nominatim_cache(i,2), Core_Utils.getCountryISO3166(nominatim_cache(i,3)));
            end
            if nargin == 1
                try
                    fh = fopen(file_path, 'wt');
                    fwrite(fh, str);
                    fclose(fh);
                catch ex
                    Core_Utils.printEx(ex);
                end
            end
        end

        function [location_info] = getLocationInfo(lat, lon, type)
            % Get reverse nominatim
            %
            % INPUT (variable)
            %   lat,lon             coordinates in degrees
            %   type                mapquest (default) / openstreetmap
            %
            % OUTPUT
            %   nominatim struct
            
            if nargin == 2
                type = 'auto';
            end
            switch type
                case {'mapquest'}
                    rrm = Remote_Resource_Manager.getInstance;
                    KEY = rrm.getMapQuestAPI;
                    if isempty(KEY)
                        Core.getLogger.addWarning('To use mapquestapi you need to have a KEY in your credentials.txt')
                        % Nominatim.openstreetmap
                        fprintf('Requesting location info of (%.4f, %.4f) from openstreetmap.org\n', lat, lon);
                        request_url = sprintf(['https://nominatim.openstreetmap.org/reverse?' ...
                            'lat=%f&' ...
                            'lon=%f&' ...
                            'format=json&' ...
                            'addressdetails=1&' ...
                            'namedetails=0&' ...
                            'extratags=0'], lat, lon);
                    else
                        fprintf('Requesting location info of (%.4f, %.4f) from mapquestapi.com\n', lat, lon);
                        request_url = sprintf(['http://open.mapquestapi.com/nominatim/v1/reverse.php?' ...
                            'key=%s&', ...
                            'format=json&' ...
                            'lat=%f&' ...
                            'lon=%f&'], KEY, lat, lon);
                    end
                case {'openstreetmap', 'auto'}
                    % Nominatim.openstreetmap
                    request_url = sprintf(['https://nominatim.openstreetmap.org/reverse?' ...
                        'lat=%f&' ...
                        'lon=%f&' ...
                        'format=json&' ...
                        'addressdetails=1&' ...
                        'namedetails=0&' ...
                        'extratags=0'], lat, lon);
            end
            try
                % First approach but somehow sometime it display log on screen
                %options = weboptions; options.CertificateFilename=(''); location_info = webread(request_url,options, 'Timeout',5);
                % Second approach:
                location_info = jsondecode(char(Core_Utils.urlread_auth(request_url,'anonymous','')));
                % if this test fail (raise an exception) it meas that there no coordinates have been found!
                test = length(location_info.address.country_code) == 2;
            catch ex
                if strcmp(type, 'auto')
                    try
                        [location_info] = Core_Utils.getLocationInfo(lat, lon, 'mapquest');
                    catch ex
                        location_info = [];
                    end
                end
            end
        end
        
        function [s, info] = urlread_auth(url, user, password)
            %URLREAD_AUTH Like URLREAD, with basic authentication
            %
            % [s,info] = urlread_auth(url, user, password)
            %
            % Returns bytes. Convert to char if you're retrieving text.
            %
            % Examples:
            % sampleUrl = 'http://browserspy.dk/password-ok.php';
            % [s,info] = urlread_auth(sampleUrl, 'test', 'test');
            % txt = char(s)
            % Matlab's urlread() doesn't do HTTP Request params, so work directly with Java
            %
            % Code found here: https://it.mathworks.com/matlabcentral/answers/376318-webread-certificate-problems-roll-your-own-java-works-webread-doesn-t-why
            % Author's unknown  
            jUrl = java.net.URL(url);
            conn = jUrl.openConnection();
            conn.setRequestProperty('Authorization', ['Basic ' matlab.net.base64encode([user ':' password])]);
            conn.connect();
            info.status = conn.getResponseCode();
            info.errMsg = char(readstream(conn.getErrorStream()));
            s = readstream(conn.getInputStream());
            
            function out = readstream(inStream)
                %READSTREAM Read all bytes from stream to uint8
                try
                    import com.mathworks.mlwidgets.io.InterruptibleStreamCopier;
                    byteStream = java.io.ByteArrayOutputStream();
                    isc = InterruptibleStreamCopier.getInterruptibleStreamCopier();
                    isc.copyStream(inStream, byteStream);
                    inStream.close();
                    byteStream.close();
                    out = typecast(byteStream.toByteArray', 'uint8'); %'
                catch err
                    out = []; %HACK: quash
                end
            end
        end
        
        
        function idx = findMO(find_list, to_find_el)
            % find the postion of the elements of to_find_el into find_list
            % find list should have unique elements
            idx = zeros(size(to_find_el));
            for i = 1: length(to_find_el)
                idx(i) = find(find_list == to_find_el(i),1);
            end
        end
        
        function [amb_idx, n_amb] = remEmptyAmbIdx(amb_idx, n_amb)
            % remove emtpy amb_idx
            %
            % SYNTAX:
            % amb_idx = Core_Utils.remEmptyAmbIdx(amb_idx, <n_amb>)
            if nargin < 2
                n_amb = max(max(amb_idx));
            end
            i = 1;
            while (i <= n_amb)
                n_ep_amb = sum(sum(amb_idx == i));
                if n_ep_amb == 0
                    n_amb = n_amb - 1;
                    amb_idx(amb_idx > i) = amb_idx(amb_idx > i) - 1;
                else
                    i = i + 1;
                end
            end
        end
        
        function num = round_even(num)
            num = round((num-2)/2)*2+2;
        end
        
        function num = round_odd(num)
            num = round((num-1)/2)*2+1;
        end
        
        function r = xcorr(x)
            % compute cross correlation
            %
            % SYNTAX:
            % xcorr = Core_Utils.xcorr(x)
            %
            % NOTE:
            % thank you Amro https://stackoverflow.com/questions/3949324/calculate-autocorrelation-using-fft-in-matlab
            len = length(x);
            
            % autocorrelation
            nfft = 2^nextpow2(2*len-1);
            r = ifft( fft(x,nfft) .* conj(fft(x,nfft)) );
            % rearrange and keep values corresponding to lags: -(len-1):+(len-1)
            r = [r(end-len+2:end) ; r(1:len)];
        end
        
        function [s,n] = getSemivariogram1D(x,mode)
            % compute 1 d semivariogram
            %
            % SYNTAX:
            %     s = Core_Utils.getSemivariogram1D(x)
            if nargin < 2
                mode = 'mean';
            end
            max_lag = length(x)-1;
            s = nan(max_lag,1);
            n = zeros(max_lag,1);
            if strcmpi(mode,'mean')
                for l = 1 : max_lag
                    diffs = (x((l+1):end) - x(1:(end-l))).^2;
                    s(l) = mean(diffs,'omitnan')/2;
                    n(l) = sum(~isnan(diffs));
                end
            elseif strcmpi(mode,'fft')
                L = length(x);
                %                 Y = fft(x)
                %                 Ae = abs(Y);
                A = Core_Utils.getSpectrum(x);
                cova = ifft(([A; flipud(A(2:end))]).^2);
                s = var(x) - cova(1:101);
            else
                for l = 1 : max_lag
                    s(l) = median((x((l+1):end) - x(1:(end-l))).^2,'omitnan')/2;
                end
            end
        end
        
        function [s,n] = getSemivariogram1DMF(x,mf)
            % compute 1 d semivariogram
            %
            % SYNTAX:
            %     s = Core_Utils.getSemivariogram1D(x)
            
            max_lag = length(x)-1;
            s = nan(max_lag,1);
            n = zeros(max_lag,1);
            
            for l = 1 : max_lag
                diffs = ((x((l+1):end) - x(1:(end-l))) ./((mf((l+1):end) + mf(1:(end-l))))*2).^2;
                s(l) = mean(diffs,'omitnan')/2;
                n(l) = sum(~isnan(diffs));
            end
            
        end
        
        function [s,n] = getSemivariogram1DMFMV(x,mf,lat,lon);
            % compute 1 d semivariogram
            %
            % SYNTAX:
            %     s = Core_Utils.getSemivariogram1D(x)
            
            max_lag = length(x)-1;
            s = nan(max_lag,1);
            n = zeros(max_lag,1);
            ang_dist = @(lat1,lat2,lon1,lon2)  acos(sin(lat1).*sin(lat2)+cos(lat1).*cos(lat2).*cos(lon1-lon2));
            for l = 1 : max_lag
                ads = ang_dist(lat((l+1):end),lat(1:(end-l)) ,lon((l+1):end),lon(1:(end-l)));
                ads = (ads/min(ads)).^(5/3);
                diffs = ((x((l+1):end) - x(1:(end-l))) ./((mf((l+1):end) + mf(1:(end-l))))*2./ads).^2;
                s(l) = mean(diffs,'omitnan')/2;
                n(l) = sum(~isnan(diffs));
            end
            
        end
        
        function [s,n] = getSemivariogram1DCS(x,cs_lid)
            % compute 1 d semivariogram accoutnign for cycle slips
            %
            % SYNTAX:
            %     s = Core_Utils.getSemivariogram1D(x)
            cs = unique([1; find(cs_lid); length(x)]);
            max_lag = length(x)-1;
            
            s = nan(max_lag,1);
            n = zeros(max_lag,1);
            for c = 2 : length(cs)
                if sum(~isnan(x(cs(c-1):cs(c)))) > 1
                    [st,nt] = Core_Utils.getSemivariogram1D(x(cs(c-1):(cs(c)-1)));
                    ntt = n(1:length(st))+nt;
                    s(1:length(st)) = zero2nan( ( nan2zero(s(1:length(st))).*n(1:length(st)) + nan2zero(st).*nt )./ntt );
                    n(1:length(st)) = ntt;
                end
            end
        end
        
        function [a,b] = logLogLineEst(y,lims)
            % compute slope and intercept in log log plane
            %
            % SYNTAX:
            %     [a,b] = Core_Utils.logLogLineEst(y,lims)
            x = (1:length(y))';
            x = x(lims(1):lims(2));
            y = y(lims(1):lims(2));
            A = [log(x) ones(size(x))];
            y = log(y);
            est = A\y;
            a=est(1);
            b=est(2);
        end
        
        function [gh11,nh11] = variofl(x1,icode)
            % SOURCE : Marcotte, Denis. "Fast variogram computation with FFT." Computers & Geosciences 22.10 (1996): 1175-1186.
            %
            % function [gh11,nh11l=variofl(x1,icode);
            %
            % function to compute variograms or covariograms, in 1D or 2D
            % the data are on a (possibly incomplete) regular grid.
            % the program computes variograms in the frequency domain by
            % using 2D-FFT.
            %
            % input: x1: data matrix. Missing values are indicated by NaN
            %
            %icode: a code to indicate which function to compute
            %=l : variogram
            %
            % =2 : covariogram
            %
            %
            % gh11: variogram or covariogram depending on icode.
            % output:
            % nh11: number of pairs available
            %
            %
            % this program uses the functions FFT2, IFFTZ, FFTlSHIFT and CONJ which are
            % standard MATLAB functions.
            [n,p]=size(x1);
            nrows=2*n-1;
            ncols=2*p-1;
            % dimensions of data matrix
            % find the closest multiple of 8 to obtain a good compromise between
            % speed (a power of 2) and memory required
            approx = 8;
            nr2=ceil(nrows/approx)*approx;
            nc2=ceil(ncols/approx)*approx;
            % form an indicator matrix:
            %      l's for all data values
            %      O's for missing values
            %
            % in data matrix, replace missing values by 0;
            x1id=~isnan(x1);
            x1(~x1id)=zeros(sum(sum(-x1id)),1); % 1 for a data value; 0 for missing
            % missing replaced by 0
            fx1=fft2(x1,nr2,nc2); % fourier transform of x1
            if icode==1
                fx1_x1=fft2(x1.*x1,nr2,nc2);
            end
            clear x1;
            fx1id=fft2(x1id,nr2,nc2);
            clear x1id
            % fourier transform of x1*x1
            % fourier transform of the indicator matrix
            % compute number of pairs at all lags
            nh11=round(real(ifft2(conj(fx1id).*fx1id)));
            % compute the different structural functions according to icode
            if icode==1
                % variogram is computed
                gh11=real(ifft2(conj(fx1id).*fx1_x1+conj(fx1_x1).*fx1id-2*conj(fx1).*fx1));
                gh11=gh11./max(nh11,1)/2;
            else
                % covariogram is computed
                ml=real(ifft2(conj(fx1).*fx1id))./max(nh11,1);
                m2=real(ifft2(conj(fx1id).*fx1))./max(nh11,1);
                clear fx1id
                gh11=real(ifft2(conj(fx1).*fx1));
                gh11=gh11./max(nh11,1)-ml.*m2;
            end
            % compute tail mean
            % compute head mean
            clear fx1 fx1id fx1_fx1
            % reduce matrix to required size and shift so that the 0 lag appears at the center of each matrix
            nh11 = [nh11(1:n,1:p) nh11(1:n,nc2-p+2:nc2); nh11(nr2-n+2,1:p) nh11(nr2-n+2:nr2,nc2-p+2:nc2)];
            gh11 = [gh11(1:n,1:p) gh11(1:n,nc2-p+2:nc2); gh11(nr2-n+2,1:p) gh11(nr2-n+2:nr2,nc2-p+2:nc2)];
            
            %             gh11=fftshift(gh11);
            %             nh11=fftshift(nh11);
        end
        
        function s = remPeriod(s,p)
            % remove signal for a given period
            %
            % SYNTAX:
            %     s = Core_Utils.remPeriod(s,p)
            ph = (1:numel(s))'/p*2*pi;
            A = [cos(ph) sin(ph)];
            N = A'*A; % less stable but faster
            B = A'*s;
            x = N\B;
            s = s - A*x;
        end
        
        function num = code2Char2Num(str2)
            % Convert a 2 char string into a numeric value (float)
            % SYNTAX
            %   num = Core_Utils.code3ch2Num(str3);
            
            num = str2(:,1:2) * [2^8 1]';
        end
        
        function str2 = num2Code2Char(num)
            % Convert a numeric value (float) of a 2 char string
            % SYNTAX
            %   str2 = Core_Utils.num2Code2ch(num)
            num = uint16(num);
            str2 = char(zeros(numel(num), 2));
            str2(:,1) = char(bitand(num, uint16(sum(2.^(8 + (0:7))))) / 2^8);
            str2(:,2) = char(bitand(num, uint16(sum(2.^(0:7)))));
        end
        
        function num = code3Char2Num(str3)
            % Convert a 3 char string into a numeric value (float)
            % SYNTAX
            %   num = Core_Utils.code3ch2Num(str3);
            
            num = str3(:,1:3) * [2^16 2^8 1]';
        end
        
        function str3 = num2Code3Char(num)
            % Convert a numeric value (float) of a 3 char string
            % SYNTAX
            %   str3 = Core_Utils.num2Code3ch(num)
            num = uint32(num);
            str3 = char(zeros(numel(num), 3));
            str3(:,1) = char(bitand(num, uint32(sum(2.^(16 + (0:7))))) / 2^16);
            str3(:,2) = char(bitand(num, uint32(sum(2.^(8 + (0:7))))) / 2^8);
            str3(:,3) = char(bitand(num, uint32(sum(2.^(0:7)))));
        end
        
        function num = code4Char2Num(str4)
            % Convert a 4 char string into a numeric value (float)
            % SYNTAX
            %   num = Core_Utils.code4ch2Num(str4);            
            num = zeros(size(str4,1),1, 'uint32');
            for i = 1:size(str4,1)
                num(i) = sum(uint64(str4(i,1:4)') .* uint64([2^24 2^16 2^8 1])');
            end
        end
        
        function num = code7Char2Num(str7)
            % Convert a 7 char string into a numeric value (float)
            % SYNTAX
            %   num = Core_Utils.code7ch2Num(str7);
            num = zeros(size(str7,1),1, 'uint64');
            for i = 1:size(str7,1)
                num(i) = sum(uint64(str7(i,1:7)') .* uint64([2^48 2^40 2^32 2^24 2^16 2^8 1])');
            end
        end

        function str4 = num2Code4Char(num)
            % Convert a numeric value (float) of a 4 char string
            % SYNTAX
            %   str4 = Core_Utils.num2Code4Char(num)
            num = uint32(num);
            str4 = char(zeros(numel(num), 4));
            str4(:,1) = char(bitand(num, uint32(sum(2.^(24 + (0:7))))) / 2^24);
            str4(:,2) = char(bitand(num, uint32(sum(2.^(16 + (0:7))))) / 2^16);
            str4(:,3) = char(bitand(num, uint32(sum(2.^(8 + (0:7))))) / 2^8);
            str4(:,4) = char(bitand(num, uint32(sum(2.^(0:7)))));
        end

        function str7 = num2Code7Char(num)
            % Convert a numeric value (float) of a 7 char string
            % SYNTAX
            %   str7 = Core_Utils.num2Code7Char(num)
            num = uint64(num);
            str7 = char(zeros(numel(num), 4));
            str7(:,1) = char(bitand(num, uint64(sum(2.^(48 + (0:7))))) / 2^48);
            str7(:,2) = char(bitand(num, uint64(sum(2.^(40 + (0:7))))) / 2^40);
            str7(:,3) = char(bitand(num, uint64(sum(2.^(32 + (0:7))))) / 2^32);
            str7(:,4) = char(bitand(num, uint64(sum(2.^(24 + (0:7))))) / 2^24);
            str7(:,5) = char(bitand(num, uint64(sum(2.^(16 + (0:7))))) / 2^16);
            str7(:,6) = char(bitand(num, uint64(sum(2.^(8 + (0:7))))) / 2^8);
            str7(:,7) = char(bitand(num, uint64(sum(2.^(0:7)))));
        end

        function str4 = unique4ch(str4)
            % Perform unique on an array of 4 char codes
            %
            % SYNTAX
            %   str4 = Core_Utilis.unique4ch(str4)
            str4 = Core_Utils.num2Code4Char(unique(Core_Utils.code4Char2Num(str4)));
        end
        
        function str3 = unique3ch(str3)
            % Perform unique on an array of 3 char codes
            %
            % SYNTAX
            %   str3 = Core_Utilis.unique3ch(str3)
            str3 = Core_Utils.num2Code3Char(unique(Core_Utils.code3Char2Num(str3)));
        end
        
        function str2 = unique2ch(str2)
            % Perform unique on an array of 2 char codes
            %
            % SYNTAX
            %   str2 = Core_Utilis.unique2ch(str3)
            str2 = Core_Utils.num2Code2Char(unique(Core_Utils.code2Char2Num(str2)));
        end
        
        function [f_status_lst, aria_err_code] = aria2cDownloadUncompress(file_name_lst, f_ext_lst, f_status_lst, date_list, out_dir)
            % Try to download files using aria2C
            %
            % INPUT
            %   file_name_list      list of file_names to download (remote path)   [cell]
            %   f_ext_lst           extension of compression ('' is valid)         [cell]
            %   f_status_lst        bool array of files existing                   [bool]
            %   date_list           GPS_Time of days of interest                   [GPS_Time]
            %   out_dir             path to out folder                             [char]
            %
            % SYNTAX
            %   f_status_lst = Core_Utils.aria2cDownloadUncompress(file_name_lst, f_ext_lst, f_status_lst, <date_list>, <out_dir>)
            %
            
            aria_err_code = 0;
            log = Core.getLogger();
            fnp = File_Name_Processor();
            rm = Remote_Resource_Manager.getInstance;
            state = Core.getCurrentSettings();
            if ispc()
                aria2c_path = '.\utility\thirdParty\aria2-extra\aria2_win\aria2c.exe';
            elseif ismac()
                aria2c_path = '/usr/local/bin/aria2c';
                if ~exist(aria2c_path, 'file')
                    aria2c_path = '/opt/homebrew/bin/aria2c';
                end
            else % is linux
                aria2c_path = '/usr/local/bin/aria2c';
                if ~exist(aria2c_path, 'file')
                    aria2c_path = '/usr/bin/aria2c';
                end
            end
            
            if ~exist(aria2c_path, 'file')
                Core.getLogger.addError('aria2c is not working, is it installed?');
                f_status_lst = false(size(file_name_lst));
                aria_err_code = 1;
                return
            end

            credentials = '';

            % Get file list path for all the files present in the list
            idf = find(~f_status_lst);
            fnl = file_name_lst(idf);
            fel = f_ext_lst(idf); % file extension list
            odl = {}; % out_dir_list
            ffp = {}; % final file path
            for i = 1 : numel(fnl)
                file_name = fnl{i};
                server = regexp(file_name,'(?<=\?{)\w*(?=})','match', 'once'); % saerch for ?{server_name} in paths
                if ~isempty(server)
                    file_name = strrep(file_name,['?{' server '}'],'');
                    [s_ip, port, user,passwd] = rm.getServerIp(server);
                    switch port
                        case {'21'} 
                            fnl{i} = ['ftp://' s_ip ':' port file_name fel{i}];
                            if ~isempty(user) && ~isempty(passwd)
                                credentials = sprintf('--ftp-user=%s  --ftp-passw=%s',user,passwd);
                            end
                        case {'443'}
                            fnl{i} = ['https://' s_ip file_name fel{i}];
                            if ~isempty(user) && ~isempty(passwd)
                                credentials = sprintf('--http-user=%s  --http-passw=%s',user,passwd);
                            end
                        otherwise
                            fnl{i} = ['http://' s_ip ':' port file_name fel{i}];
                            if ~isempty(user) && ~isempty(passwd)
                                credentials = sprintf('--http-user=%s  --http-passw=%s',user,passwd);
                            end
                    end
                else
                    fnl{i} = [file_name fel{i}];
                end
                if nargin < 5
                    out_dir = Core.getState.getFileDir(file_name);
                end
                if nargin >= 4 && ~isempty(date_list)
                    out_dir = fnp.dateKeyRep(out_dir, date_list.getEpoch(date_list.length - idf(i) + 1),'',state.vmf_res,upper(state.vmf_source));
                end
                odl{i} = out_dir;
                [~, name, ext] = fileparts(fnl{i});
                if strcmp(ext,'.Z') || strcmp(ext,'.gz')
                    ffp{i} = fullfile(out_dir, name);
                else
                    ffp{i} = fullfile(out_dir, [name ext]);
                end
            end
            
            % if I have at least one file to download
            if numel(odl) > 0
                i = 0;
                file_name = fullfile('.', 'reserved', sprintf('tmpAriaDownload_%s.lst', datestr(now, 'yyyymmddHHMMSSFFF')));
                if (exist(['.' filesep 'reserved' filesep], 'dir') == 0)
                    mkdir(['.' filesep 'reserved']);
                end
                fid = fopen(file_name, 'Wb');
                if fid < 0
                    log.addWarning(['Writing on "' file_name '" is not possible' char(10) 'aria2 could not work']);
                else
                    str = '';
                    old_od = odl{1};
                    while i <= numel(fnl)
                        i = i + 1;
                        if (i < numel(fnl)) && strcmp(odl{i}, old_od)
                            str = sprintf('%s%s\n', str, fnl{i});
                        else
                            
                            if i <= numel(fnl)
                                if i == numel(fnl) && strcmp(odl{i}, old_od)
                                    str = sprintf('%s%s\n', str, fnl{i});
                                end
                                
                                % call aria
                                % check for .Z or .gz compressed files too
                                fwrite(fid, str, '*char');
                                fclose(fid);
                                if ~isempty(str)
                                    if ~exist(old_od, 'file')
                                        mkdir(old_od);
                                    end
                                    log.addMessage(sprintf('Executing \n  aria2c -c -i %s -d %s\n  File download list:', file_name, old_od));
                                    log.addMessage(log.indent(sprintf('%s', str)));
                                    try
                                       if any(strfind(str, 'ftp://garner')) || any(strfind(str, 'ftp://nfs.kasi'))
                                            str_parallel = '1';
                                        elseif any(strfind(str, 'ftp://gssc.esa')) 
                                            str_parallel = '2';
                                        else
                                            str_parallel = '20';
                                        end
                                        if ispc()
                                            %dos(sprintf('"%s" %s -j %s -c -i %s -d %s >nul 2>&1', aria2c_path, credentials, str_parallel, file_name, old_od)); % suppress output
                                            dos(sprintf('"%s" %s -j 1 -c -i %s -d %s', aria2c_path, credentials, file_name, old_od)); % do not suppress output
                                        else
                                            %dos(sprintf('%s %s -j %s -c -i %s -d %s &> /dev/null', aria2c_path, credentials, str_parallel, file_name, old_od));  % suppress output
                                            dos(sprintf('%s %s -j 1 -c -i %s -d %s ', aria2c_path, credentials, file_name, old_od));  % do not suppress output
                                        end
                                    catch
                                        f_status_lst = false(size(file_name_lst));
                                        aria_err_code = 1;
                                        log.addError('aria2c is not working, is it installed?');
                                    end
                                end
                                % Check for zero byte files (download errors)
                                for f = 1 : numel(fnl)
                                    [~, out_file_name, out_file_ext] = fileparts(fnl{f});
                                    out_file_path = [old_od, filesep, out_file_name, out_file_ext];
                                    if exist(out_file_path, 'file') == 2
                                        file_info = dir(out_file_path);
                                        if file_info.bytes == 0
                                            % the file is empty
                                            delete(out_file_path)
                                            log.addError(sprintf('%s download failed\nThe file is probably missing', [out_file_name, out_file_ext]));
                                        end
                                    end
                                end
                                % open file list for the next set
                                fid = fopen(file_name, 'Wb');
                                str = '';
                            end
                        end
                        if i <= numel(fnl)
                            decrement = 0;
                            if ~(strcmp(odl{i}, old_od))
                                decrement = 1; % this last file have not been downloaded!
                            end
                            old_od = odl{i};
                            i = i - decrement;
                        end
                    end
                    fclose(fid);
                    delete(file_name);
                end
                
                % once the file have been downloaded decompress and test presence
                for i = 1 : numel(fnl)
                    if strcmp(fel{i},'.Z') || strcmp(fel{i},'.gz')
                        if (isunix())
                            system(['gzip -d -f "' ffp{i} fel{i} '" &> /dev/null']);
                        else
                            try                                        
                                [app_path] = Core.getInstallDir;
                                [status, result] = system(['"' app_path '\utility\thirdParty\7z1602-extra\7za.exe" -y x "'   ffp{i} fel{i} '" -o"'  odl{i} '"']); 
                                if (status == 0)
                                    status = true;
                                elseif (status > 1)
                                    this.log.addError(sprintf('Please decompress the %s file before trying to use it in Breva!!!\nError 7za: %s', [ffp{i} fel{i}], result));
                                end
                                delete([ffp{i} fel{i}]);
                            catch
                                this.log.addError(sprintf('Please decompress the %s file before trying to use it in Breva!!!', [ffp{i} fel{i}]));
                                status = false;
                            end
                        end
                    end
                end
                
                % update f_status_lst
                for i = 1 : numel(ffp)
                    f_status_lst(idf(i)) = exist(ffp{i}, 'file') == 2;
                end
            end
            
        end
        
        function [status] = downloadHttpTxtResUncompress(filename, out_dir, user, passwd, https_flag)
            if nargin < 3
                user = '';
                passwd = '';
            end
            if nargin < 5
                https_flag = false;
            end
            
            log = Core.getLogger();
            fnp = File_Name_Processor();
            try


                auth_str = [user ':' passwd];
                base64Str = ['Basic ' matlab.net.base64encode(auth_str)];
                headers = {'Authorization', base64Str};
                options = weboptions('HeaderFields',headers);
                options.ContentType = 'text';
                options.Timeout = 5;
                options.Username = user;
                options.Password = passwd;
                [remote_location, filename, ext] = fileparts(filename);
                filename = [filename ext];
                log.addMessage(log.indent(sprintf('downloading %s ...',filename)));
                compressed_name = '';
                status = true;
                if ~isempty(out_dir) && ~exist(out_dir, 'dir')
                    mkdir(out_dir);
                end
                   
                if any(strfind(remote_location, 'cddis'))
                    % cddis download data even if the compression estension is not specified, but with the wrong extension name!
                    % I need to force the extension
                    if contains(remote_location, 'igs20') || length(filename) > 20
                        filename = [filename '.gz'];
                    else
                        filename = [filename '.Z'];
                    end
                    compressed_name = filename;
                end
                if strcmp(ext, '.zip')
                    compressed_name = filename;
                end
                
                try
                    if ~https_flag
                        txt = websave(fullfile(out_dir, filename), ['http://' remote_location '/' filename], options);
                    else
                        txt = websave(fullfile(out_dir, filename), ['https://' remote_location '/' filename], options);
                    end
                catch ex
                    if any(strfind(remote_location, 'cddis'))
                        % cddis download data even if the compression estension is not specified, but with the wrong extension name!
                        % Remove Z extension
                        filename = filename(1:end-2);
                        compressed_name = '';
                    end
                    if instr(ex.message, '404')
                        try
                            compressed_name = [filename, '.gz'];
                            if ~https_flag
                                txt = websave(fullfile(out_dir, compressed_name), ['http://' remote_location '/' compressed_name], options);
                            else
                                txt = websave(fullfile(out_dir, compressed_name), ['https://' remote_location '/' compressed_name], options);
                            end
                        catch ex
                            if instr(ex.message, '404')
                                try
                                    compressed_name = [filename, '.Z'];
                                    if ~https_flag
                                        txt = websave(fullfile(out_dir, compressed_name), ['http://' remote_location '/' compressed_name], options);
                                    else
                                        txt = websave(fullfile(out_dir, compressed_name), ['https://' remote_location '/' compressed_name], options);
                                    end
                                catch
                                    status = false;
                                end
                            end
                        end
                    elseif instr(ex.message, '401')
                        status = false;
                        log.addError('Unauthorized, please add credentials to credentials.txt');
                    end
                end
                if contains(txt,'.html')
                    status = false;
                    delete(txt);
                    log.addError(sprintf('Unauthorized, please check credentials into "%s"', Core.getFilePath('credentials', false)));
                end
                if status
                    status = false;
                    if ~isempty(compressed_name)
                        compressed_name = fnp.checkPath(fullfile(out_dir, compressed_name));
                        ok_status = Core_Utils.uncompressFile(compressed_name);
                    end
                    status = true;
                    log.addMessage(' Done');
                end
            catch ex
                status = false;
            end
        end

        function ok_status = uncompressFile(full_file_path)
            % Uncompress inplace a compressed file
            % 
            % SYNTAX
            %   ok_status = Core_Utils.uncompressFile(full_file_path);
            ok_status = true;
            [out_dir, ~, fext] = fileparts(full_file_path);
            if strcmp(fext,'.Z') || strcmp(fext,'.gz')  || strcmp(fext,'.zip')
                if (isunix())
                    if strcmp(fext,'.Z') || strcmp(fext,'.gz')
                        system(['gzip -d -f "' full_file_path '" &> /dev/null &']);
                    else
                        res = system(['unzip "' full_file_path '" -d "' out_dir '" &> /dev/null &']);
                        if res == 0
                            % delete zip file
                            delete(full_file_path);
                        end
                    end
                else
                    try
                        [app_dir] = Core.getInstallDir();
                        [ok_status, result] = system(['"' app_dir '\utility\thirdParty\7z1602-extra\7za.exe" -y x "' full_file_path '" -o"'  out_dir '"']); %#ok<ASGLU>
                        if (ok_status == 0)
                            ok_status = true;
                            delete(full_file_path);
                        end
                    catch
                        this.log.addError(sprintf('Please decompress the %s file before trying to use it in Breva!!!', full_file_path));
                        ok_status = false;
                    end
                end
            end
        end

        function ok_status = compressFile(full_file_path)
            % Compress inplace an uncompressed file
            % Uses gz compression
            %
            % SYNTAX
            %   ok_status = Core_Utils.compressFile(full_file_path);
            ok_status = true;
            log = Logger.getInstance;
            if exist([full_file_path '.gz'], 'file')
                delete([full_file_path '.gz'])
            end
            if (isunix())
                system(['gzip "' full_file_path '" &> /dev/null &']);
            else
                log = Logger.getInstance;
                try
                    [app_path] = Core.getInstallDir();
                    [app_dir] = fileparts(app_path);
                    [ok_status, result] = system(['"' app_dir '\utility\thirdParty\7z1602-extra\7za.exe" a "' full_file_path '.gz" "' full_file_path '"']); %#ok<ASGLU>
                    if (ok_status == 0)
                        ok_status = true;
                        log.addStatusOk(sprintf('File "%s" compressed successfully', full_file_path));
                    else
                        log.addStatusOk(sprintf('File "%s" compressed successfully', full_file_path));
                    end

                catch
                    log.addError(sprintf('gzip compression of %s is not working', full_file_path));
                    ok_status = false;
                end
            end
            if (ok_status)
                log.addStatusOk(sprintf('File "%s" compressed successfully', full_file_path));
            else
                log.addError(sprintf('File "%s" was not compressed', full_file_path));
            end
        end
        
        function [status, ext] = checkHttpTxtRes(filename, user, passwd)
            ext = '';
            if isunix()
                if any(strfind(filename, 'cddis'))
                    status = true; % To check cddis you need credentials, so here it is not implemented
                else
                    if ismac()
                        % Sometimes mac users does not have wget'
                        rem_check_cmd = 'curl --head ';
                        if nargin > 1 && ~isempty(user)
                            rem_check_cmd = [rem_check_cmd sprintf('---insecure -u %s:%s ', user, passwd)];
                        end
                    else
                        % curl seems to have some problems with matlb libraries
                        % under some Linux installations, switching to wget --spyder
                        rem_check_cmd = 'wget --spider --no-hsts --no-check-certificate ';
                        if nargin > 1 && ~isempty(user)
                            rem_check_cmd = [rem_check_cmd sprintf('--user=%s --password=%s ', user, passwd)];
                        end
                    end
                    
                    [resp, txt] = system([rem_check_cmd filename]);
                    if ~isempty(strfind(txt,' 200 OK')) || ~isempty(strfind(txt,' 302 Found')) || ~isempty(strfind(txt,' 302 Moved Temporarily')) %#ok<STREMP>
                        status = true;
                    else
                        [resp, txt] = system([rem_check_cmd filename '.gz']);
                        if ~isempty(strfind(txt,' 200 OK')) || ~isempty(strfind(txt,' 302 Found')) || ~isempty(strfind(txt,' 302 Moved Temporarily')) %#ok<STREMP>
                            ext = '.gz';
                            status = true;
                        else
                            [resp, txt] = system([rem_check_cmd filename '.Z']);
                            if ~isempty(strfind(txt,' 200 OK')) || ~isempty(strfind(txt,' 302 Found')) || ~isempty(strfind(txt,' 302 Moved Temporarily')) %#ok<STREMP>
                                ext = '.Z';
                                status = true;
                            else
                                status = false;
                            end
                        end
                    end
                end
            else
                log = Logger.getInstance;
                log.addWarning('HTTP check is implemented only for Unix systems')
                status = true; % !!! to be implemented
            end
        end
        
        function station_list = getStationList(dir_path_list, file_ext, flag_recursive)
            % Get the list of stations present in a folder (with keys substituted)
            %
            % SYNTAX
            %   station_list = Core_Utilis.getStationList(dir_path)
            
            if nargin == 1
                file_ext = '.';
            else
                file_ext = ['[' file_ext ']'];
            end
            
            if nargin == 3 && flag_recursive
                dir_path_list = strsplit(genpath(dir_path_list), ':');
            else
                dir_path_list = {dir_path_list};
            end
            
            recursive_dirs = {};
            for d = 1 : numel(dir_path_list)
                dir_path = dir_path_list{d};
                if ~isempty(dir_path)
                    try
                        % Calling dos is faster than dir with large directories
                        if isunix
                            [~, d] = dos(['ls "' dir_path '"']);
                            dir_list = strsplit(d);
                            dir_list = dir_list(1:end-1);
                        else
                            [~, d] = dos(['dir "' dir_path '"']);
                            dir_list = strsplit(d);
                            dir_list = dir_list(1:end-1);
                        end
                    catch
                        dir_list = dir(dir_path);
                        dir_list = {dir_list.name};
                    end
                    for i = 1 : numel(dir_list)
                        dir_list{i} = [dir_path(length(dir_path_list{1})+2:end) filesep dir_list{i}];
                    end
                end
                recursive_dirs = [recursive_dirs dir_list];
            end
            dir_list = recursive_dirs;
            
            % search for station files STAT${DOY}${S}${QQ}.${YY}
            file_list = {};
            for d = 1 : numel(dir_list)
                tmp = strsplit(dir_list{d}, filesep);
                file_name = tmp{end};
                file_name_len = numel(file_name);
                rin2_start = regexp(dir_list{d}, ['.{4}[0-9]{3}.{1}[0-9]{2}[\.]{1}[0-9]{2}' file_ext '{1}'], 'once');
                rin3_start = regexp(dir_list{d}, '\_[0-9]{4}[0-9]{3}[0-9]{4}\_', 'once');
                if (file_name_len == 14) && ~isempty(rin2_start)
                    file_list = [file_list; {[dir_list{d}(1:rin2_start) '${DOY}${S}${QQ}.${YY}' dir_list{d}(end)]}];
                elseif ~isempty(rin3_start)
                    file_list = [file_list; {[dir_list{d}(1:rin3_start) '${YYYY}${DOY}' dir_list{d}(rin3_start + 8 : end)]}]; %#ok<AGROW>
                end
            end
            
            % search for station files STAT${DOY}${S}.${YY}
            for d = 1 : numel(dir_list)
                tmp = strsplit(dir_list{d}, filesep);
                file_name = tmp{end};
                file_name_len = numel(file_name);
                rin2_start = regexp(dir_list{d}, ['[0-9]{3}.{1}[\.]{1}[0-9]{2}' file_ext '{1}'], 'once');
                if (file_name_len == 12) && ~isempty(rin2_start)
                    file_list = [file_list; {[dir_list{d}(1:rin2_start-1) '${DOY}${S}.${YY}' dir_list{d}(end)]}]; %#ok<AGROW>
                    %file_list = [file_list; dir_list{d}(1:4)];
                end
            end
            station_list = unique(file_list);
        end
        
        function data = injectData(data1, data2, idx1, idx2, data_size)
            % isert data2 into data1 at the position definied by idx1 and idx2
            % idx1 - 1 is the last element of data1 to be putted before data2 (0  if none)
            % idx2 + 1 is the first element of data1 to be put after data2 (data1 length if none)
            %
            % SYNTAX
            %   data = Core_Utils.injectData(data1, data2, idx1, idx2)
            if nargin == 5
                % check data size
                [m, n] = size(data2);
                if not(m == data_size(1) && n == data_size(2))
                    data2 = nan(data_size(1), data_size(2));
                end
            end
            if ~isempty(data2)
                if isempty(data1)
                    data1 = nan(idx1 - 1, size(data2,2));
                end
                if size(data1, 1) < idx2
                    data1((end + 1) : idx2, :) = nan;
                end
                data = [data1(1 : idx1 - 1, :); data2; data1(idx2 + 1 : end, :)];
            else
                data = data1;
            end
        end
        
        function data = injectSmtData(data_lft, data_rgt, idx_smt1, idx_smt2, time1, time2, id_stop, id_start, interpolate)
            % inject smoothed data
            %
            % INPUT:
            %   data_lft : all data left
            %   data_right : all data right
            %   idx_smt1 : which data of data_lft are to be smoothed
            %   idx_smt2 : which data of data_rgt are to be smoothed
            %   time_1: time of the left data to be smoothed
            %   time_2: time of the right data to be smoothed
            %   id_stop: first epoch of time_1 that should not be kept
            %   id_start: first epoch of time_2 to keep
            %
            % SYNTAX:
            %   data = Core_Utils.injectSmtData(data_lft, data_rgt, idx_smt1, idx_smt2, time_1, time_2, id_start)
            if nargin < 9
                interpolate = true;
            end
            
            % Get the part of data to interp
            data_tosmt_lft = data_lft(idx_smt1);
            data_tosmt_rgt = data_rgt(idx_smt2);
            
            % we use mat time, is easier and we do not need extreme precision
            time1 = time1.getMatlabTime();
            time2 = time2.getMatlabTime();
            
            [idx1, idx2, time_tot] = Core_Utils.unionOrderedDouble(time1, time2, median([diff(time1); diff(time2)])/4); % 1/4 the rate tolerance
            
            mix_len = min(0.007, abs((time2(1) - time1(end)))/20); % <= empirically found
            w2 = 1 ./ (1 + exp(-((time_tot - mean(time_tot)) / mix_len))); % todo: scale to ensure [0 1]
            w1 = 1 - w2;
            n_out = size(time_tot);
            data1 = nan(n_out);
            data2 = nan(n_out);
            data1(idx1) = data_tosmt_lft;
            data2(idx2) = data_tosmt_rgt;
            %id_start = idx1(id_start);
            %id_ko = ((isnan(data1) & (1 : n_out)' < id_start) |(isnan(data2) & (1 : n_out)' >= id_start)) & ~(isnan(data1) & isnan(data2)); %?? should be used -> yes beacuse time is injected deleting overlapping times
            id_keep = unique([idx1(1:(id_stop-1)); idx2(id_start:end)]);
            % Interpolate missing data
            if interpolate
                is_nan = find(isnan(data1));
                extr_lft = is_nan(time_tot(is_nan) <= min(time1));
                extr_rgh = is_nan(time_tot(is_nan) >= max(time1));
                data1(extr_lft) = data_tosmt_lft(1);
                data1(extr_rgh) = data_tosmt_lft(end);
                if any(~isnan(data_tosmt_lft)) &&  any(isnan(data1))
                    data1 = simpleFill1D(data1, isnan(data1), 'linear', time_tot);
                end
                is_nan = find(isnan(data2));
                extr_lft = is_nan(time_tot(is_nan) <= min(time2));
                extr_rgh = is_nan(time_tot(is_nan) >= max(time2));
                data2(extr_lft) = data_tosmt_rgt(1);
                data2(extr_rgh) = data_tosmt_rgt(end);
                if any(~isnan(data_tosmt_rgt)) && any(isnan(data2))
                    data2 = simpleFill1D(data2, isnan(data2), 'linear', time_tot);
                end
            else
                data1(isnan(data1)) = data2(isnan(data1));
                data2(isnan(data2)) = data1(isnan(data2));
            end
            
            % Do not Merge nans
            last_ok = find(~isnan(data1), 1, 'last');
            if isempty(last_ok)
                last_ok = 0;
            end
            w1((last_ok + 1) : end) = 0;
            data1((last_ok + 1) : end) = 0;
            w2((last_ok + 1) : end) = 1;
            
            % Do not Merge nans
            first_ok = find(~isnan(data2), 1, 'first');
            if isempty(last_ok)
                first_ok = numel(data2) + 1;
            end
            w1(1 : (first_ok - 1)) = 1;
            w2(1 : (first_ok - 1)) = 0;
            data2(1 : (first_ok - 1)) = 0;
            
            data = w1.*data1 + w2.*data2;
            %data(id_ko) = [];
            data = data(id_keep);
            data = [data_lft(~idx_smt1); data; data_rgt(~idx_smt2)];
        end
        
        function [idx1, idx2, double_tot] = unionOrderedDouble(double_1, double_2, threshold)
            % given two ordered double give the index of the two vector in the joint vector considering the threshold
            %
            % SYNTAX
            % [idx1, idx2] = Core_Utils.unionOrderedDouble(double_1, double_2, threshold)
            l1 = length(double_1);
            l2 = length(double_2);
            idx1 = zeros(l1,1);
            idx2 = zeros(l2,1);
            i = 1;
            j = 1;
            tot = 1;
            while  i <=  l1 && j <= l2
                if abs(double_1(i) - double_2(j)) < threshold
                    idx1(i) = tot;
                    idx2(j) = tot;
                    i = i + 1;
                    j = j + 1;
                    tot = tot + 1;
                elseif double_1(i) < double_2(j)
                    idx1(i) = tot;
                    i = i + 1;
                    tot = tot + 1;
                else
                    idx2(j) = tot;
                    j = j + 1;
                    tot = tot +1;
                end
            end
            if j > l2 && i <= l1
                % idx_end = (i : l1) -l1 + tot; % WRONG
                idx_end = (i : l1) -i + tot;
                idx1(i : l1) = idx_end;
            elseif i > l1 && j <= l2
                % idx_end = (j : l2) -l2 + tot; % WRONG
                idx_end = (j : l2) -j + tot;
                idx2(j : l2) = idx_end;
            end
            double_tot = zeros(max(max(idx1), max(idx2)), 1);
            double_tot(idx1) = double_1;
            double_tot(idx2) = double_2;
        end
        
        function [ids] = findAinB(cellA,cellB)
            % find the index of cella in cellb
            %
            % SYNTAX
            %     [ids] = Core_Utils.findAinB(cellA,cellB)
            if ~iscell(cellA)
                cellA = {cellA};
            end
            lB = length(cellB);
            lA = length(cellA);
            ids = zeros(lA,1);
            for i = 1 : lA
                not_found = true;
                j = 1;
                while j <= lB && not_found
                    if strcmp(num2str(cellA{i}),num2str(cellB{j}))
                        ids(i) = j;
                        not_found = false;
                    end
                    j = j+1;
                end
            end
            %ids(ids==0) = [];
        end
        
        function [wl_cyle_out, frac_bias] = getFracBias(wl_cycle, weigth)
            % get the common frac bias between cycles
            % NOTE: very coarse/nobrain/empirical solution - > a simpler one should be found
            %
            % SYNTAX
            %   [wl_cyle, frac_bias] = Core_Utils.getFracBias(wl_cycle)
            if nargin < 2
                weigth = ones(size(wl_cycle));
            end
            frac_bias = zeros(1,3);
            wl_cycle_frac = zeros(size(wl_cycle,1),3);
            % get receiver wsb
            wl_cycle_frac(:,1) = zero2nan(wl_cycle) - floor(zero2nan(wl_cycle));
            wl_cycle_frac(:,2) = zero2nan(wl_cycle) - round(zero2nan(wl_cycle));
            wl_cycle_frac(:,3) = zero2nan(wl_cycle) - ceil(zero2nan(wl_cycle));
            frac_bias(1) = median(wl_cycle_frac(:,1),'omitnan');
            frac_bias(2) = median(wl_cycle_frac(:,2),'omitnan');
            frac_bias(3) = median(wl_cycle_frac(:,3),'omitnan');
            wl_cyle_var = mean(abs(wl_cycle_frac-repmat(frac_bias,size(wl_cycle_frac,1),1)),'omitnan');
            [~,idx] = min(wl_cyle_var);
            frac_bias = frac_bias(idx);
            wl_cycle_frac = wl_cycle_frac(:,idx);
            a = 0;
            idx_rw = ones(size(wl_cycle_frac));
            while sum(zero2nan(wl_cycle_frac - frac_bias) < -0.5) > 0  || a < 4
                wl_cycle_frac((wl_cycle_frac-frac_bias) < -0.5) = wl_cycle_frac((wl_cycle_frac-frac_bias) < -0.5) + 1;
                idx_rw(:) = 1;
                e = abs((wl_cycle_frac-frac_bias)) > 0.2;
                idx_rw(e) = 1./((wl_cycle_frac(e)-frac_bias)/0.2).^2; %idx_reweight
                idx_rw = idx_rw  .* weigth;
                idx_rw = idx_rw / sum(idx_rw);
                frac_bias = sum((wl_cycle_frac).*idx_rw);
                a = a+1;
            end
            wl_cyle_out = wl_cycle - frac_bias;
        end
        
        function [response] = timeIntersect(time1_st, time1_end, time2_st, time2_end)
            % check whether times of time bound 1 intersect with times of time bound 2
            %
            % SYNTAX
            % [response] = Core_Utils.timeIntersect(time1_st,time1_en, time2_st, time2_en)
            response = time1_st <= time2_end & time1_end >= time2_st;
            %(time2_st <= time1_st & time2_end >= time1_st) | (time2_st <= time1_en & time2_end >= time1_en);
        end
                
        function [ it, st, ilons, ilone, slon, ilat, slat] = getIntIdx(data, first_lat, dlat, first_lon, dlon, first_t, dt, lat, lon, t)
            % get  interpolating index
            [nlat , nlon, nt] = size(data);
            
            lon(lon < first_lon) = lon(lon < first_lon) + nlon * dlon; %% to account for earth circularity
            % find indexes and interpolating length
            % time
            it = max(min(floor((t - first_t)/ dt)+1,nt-1),1);
            st = max(min(t - first_t - (it-1)*dt, dt), 0) / dt;
            st = serialize(st);
            
            % lat
            ilat = max(min(floor((lat - first_lat)/ dlat)+1,nlat-1),1);
            slat = min(max(lat - first_lat - (ilat-1)*dlat, dlat), 0) / dlat;
            
            % lon
            ilons = max(min(floor((lon - first_lon)/ dlon)+1,nlon),1);
            ilone = ilons +1;
            ilone(ilone > nlon) = 1;
            slon = max(min(lon - first_lon- (ilons-1)*dlon, dlon), 0) / dlon;
        end

        function [ it, st, ilons, ilone, slon, ilat, slat] = getIntIdxRot(data, first_lat, dlat, first_lon, dlon, first_t, dt, lat, lon, t)
            % get  interpolating index
            [nlat , nlon, nt] = size(data);
            
            if isa(t, 'GPS_Time')
                t = t.getGpsTime;
            end  
            if isa(first_t, 'GPS_Time')
                first_t = first_t.getGpsTime;
            end  
            % find indexes and interpolating length
            % time
            it = max(min(floor((t - first_t)/ dt)+1,nt-1),1);
            st = max(min(t - first_t - (it-1)*dt, dt), 0) / dt;
            st = serialize(st);
            
            % lat
            ilat = max(min(floor((lat - first_lat)/ dlat)+1,nlat-1),1);
            slat = min(max(lat - first_lat - (ilat-1)*dlat, dlat), 0) / dlat;
            
            % Correct for Earth rotation 
            rot.rate = 360 * (dt/86400);    % degree per time step (dt)
            rot.prev = + rot.rate * st;     % move East (ahead in time) the past map
            rot.next = - rot.rate * (1-st);  % move West (back in time) the future map

            % Previous map lon
            lon(lon < (first_lon - rot.prev)) = lon(lon < (first_lon - rot.prev)) + nlon * dlon; %% to account for earth circularity            
            
            ilons(:,1) = max(min(floor((lon - (first_lon - rot.prev))/ dlon)+1,nlon),1);
            ilone(:,1) = ilons(:,1) + 1;
            ilone(ilone(:,1) > nlon,1) = mod(ilone(ilone(:,1) > nlon,1)+1, nlon)-1;
            slon(:,1) = max(min(lon - (first_lon - rot.prev)- (ilons(:,1)-1)*dlon, dlon), 0) / dlon;

            % Next map lon
            lon(lon < (first_lon - rot.next)) = lon(lon < (first_lon - rot.next)) + nlon * dlon; %% to account for earth circularity            
            
            ilons(:,2) = max(min(floor((lon - (first_lon - rot.next))/ dlon)+1,nlon),1);
            ilone(:,2) = ilons(:,2) + 1;
            ilone(ilone(:,2) > nlon,2) = mod(ilone(ilone(:,2) > nlon,2)+1, nlon)-1;
            slon(:,2) = max(min(lon - (first_lon - rot.next)- (ilons(:,2)-1)*dlon, dlon), 0) / dlon;
        end
        
        function response = permutedEqual(str1, str2)
            % check if the two variables are permuted version of the same sequence
            %
            % SYNTAX:
            %     response =Core_Utils.pertutedEqual(var1, var2)
            if length(str1) ~= length(str2)
                response = false;
            else
                ll = length(str1);
                found_all = true;
                i = 1;
                j = 1;
                while i <= ll
                    found = false;
                    while j <= ll && ~found
                        found = found || str1(i) == str2(j);
                        j = j+1;
                    end
                    found_all = found_all && found;
                    i = i+1;
                end
                response = found_all;
            end
            
        end
        
        function A = remBFromA(A,B)
            % reomve the lement of B from A
            for i = 1 : length(B)
                A(A==B(i)) = [];
            end
        end
        
        function [iA] = inverseByPartsDiag(A,idx1,idx2)
            % inverse by partitioning of a diagonal block matrix
            %
            % SYNTAX:
            %  [iA] = inverseByPartsDiag(A,idx1,idx2)
            sz1 = sum(idx1);
            sz2 = sum(idx2);
            A11 = A(idx1, idx1);
            A22 = A(idx2, idx2);
            A21 = A(idx2, idx1);
            A12 = A(idx1, idx2);
            iA11 = spdiags(1./diag(A11), 0, sz1, sz1);
            iA22 = spdiags(1./diag(A22), 0, sz2, sz2);
            r1 = A12*iA22*A21;
            r2 = A21*iA11*A12;
            A21 = -iA22*A21;
            A12 = -iA11*A12;
            A11 = A11 - r1;
            A22 = A22 - r2;
            iA11 = spdiags(1./diag(A11), 0, sz1, sz1);
            iA22 = spdiags(1./diag(A22), 0, sz2, sz2);
            iA = sparse(size(A));
            iA(idx1,idx1) = iA11;
            iA(idx1,idx2) = A12*iA22;
            iA(idx2,idx1) = A21*iA11;
            iA(idx2,idx2) = iA22;
        end
        
        function [iA] = inverseByRecPartsDiag(A,indices)
            % inverse by recursive partitioning of a diagonal block matrix 
            % 
            % SYNTAX:
            %    [iA] = inverseByPartsDiags(A,indices)
            n_blk = length(indices);
            
            % slit the indece in two
            part1 = 1:ceil(n_blk/2);
            part2 = (ceil(n_blk/2)+1):n_blk;
            if length(part1) > 1
                composite1 = true;
            else
                composite1 = false;
                sz1 = sum(indices{part1});
            end
            if length(part2) > 1
                composite2 = true;
            else
                composite2 = false;
                sz2 = sum(indices{part2});
            end
            
            idx1 = false(size(indices{1}));
            idx2 = false(size(indices{1}));

            for i = 1 : length(part1)
                idx1 = idx1 | indices{part1(i)};
            end
            for i = 1 : length(part2)
                idx2 = idx2 | indices{part2(i)};
            end
            
            % actual inversion
            A11 = A(idx1, idx1);
            A22 = A(idx2, idx2);
            A21 = A(idx2, idx1);
            A12 = A(idx1, idx2);
            if composite1
                new_indeces1 = {};
                for i = 1 : length(part1)
                    idx_tmp  = indices{part1(i)};
                    new_indeces1{i} = idx_tmp(idx1);
                end
                iA11 = Core_Utils.inverseByRecPartsDiag(A11,new_indeces1);
            else
                iA11 = spdiags(1./diag(A11), 0, sz1, sz1);
            end
            if composite2
                new_indeces2 = {};
                for i = 1 : length(part2)
                    idx_tmp  = indices{part2(i)};
                    new_indeces2{i} = idx_tmp(idx2);
                end
                iA22 = Core_Utils.inverseByRecPartsDiag(A22,new_indeces2);
            else
                iA22 = spdiags(1./diag(A22), 0, sz2, sz2);
            end
            r1 = A12*iA22*A21;
            r2 = A21*iA11*A12;
            A21 = -iA22*A21;
            A12 = -iA11*A12;
            A11 = A11 - r1;
            A22 = A22 - r2;
            if composite1
                iA11 = Core_Utils.inverseByRecPartsDiag(A11,new_indeces1);
            else
                iA11 = spdiags(1./diag(A11), 0, sz1, sz1);
            end
            if composite2
                iA22 = Core_Utils.inverseByRecPartsDiag(A22,new_indeces2);
            else
                iA22 = spdiags(1./diag(A22), 0, sz2, sz2);
            end
            iA = sparse(size(A));
            iA(idx1,idx1) = iA11;
            iA(idx1,idx2) = A12*iA22;
            iA(idx2,idx1) = A21*iA11;
            iA(idx2,idx2) = iA22;
            
        end
        
        function [fb, frac_b_mat]= estimateFracBias(obs_cy, cycle_slip)
            % estimate the common factional bias to all the obesravtions
            %
            % SYNTAX:
            %    fb = Network.estimateFracBias(obs_cy, cycle_slip)
            amb_idx = Core_Utils.getAmbIdx(cycle_slip, obs_cy);
            frac_cy = obs_cy;
            n_arcs = max(amb_idx(:,end));
            frac_b = zeros(n_arcs,1);
            num_ep = zeros(n_arcs,1);
            frac_b_mat = nan(size(obs_cy));
            for i = 1 : n_arcs
                idx_amb = amb_idx == i;
                amb = round(median(obs_cy(idx_amb),'omitnan'));%floor(strongMean(obs_cy(idx_amb),1, 0.95,2.5));
                frac_cy(idx_amb) = frac_cy(idx_amb) - amb;
                frac_b(i) = median(frac_cy(idx_amb),'omitnan'); %strongMean(frac_cy(idx_amb),1, 0.95,2.5);
                num_ep(i) = sum(sum(idx_amb));
                frac_b_mat(idx_amb) = frac_b(i);
            end
            frac_b_mat = frac_b_mat+0*obs_cy;
            fb = Core_Utils.circularModedRobustMean(frac_b_mat(:));
            frac_cy = nan(size(obs_cy));
            frac_b_mat = nan(size(obs_cy));
            for i = 1 : n_arcs
                idx_amb = amb_idx == i;
                amb = round(median(obs_cy(idx_amb) - fb,'omitnan'));%floor(strongMean(obs_cy(idx_amb),1, 0.95,2.5));
                frac_cy(idx_amb) = obs_cy(idx_amb) - amb - fb;
                frac_b_mat(idx_amb) = amb;
            end
            fb = fb + strongMean(frac_cy(:),1, 0.90,2.5);
        end
        
        function fr_cy = circularMean(cycle, obs_weigth)
            % estimate the mean for data over 0 -1
            %
            % SYNTAX:
            %     fr_cy = Core_Utils.circularMean(cycle)
            
            % elimnate the nan
            if nargin < 2
                obs_weigth = ones(size(cycle));
            end
            idx_nan = isnan(cycle);
            cycle(idx_nan) = [];
            obs_weigth(idx_nan) = [];
            cycle = cycle*2*pi;
            % compute the unit vectors, take the mean and recompute the
            % average
            
            unit = [cos(cycle) sin(cycle)];
            obs_weigth = obs_weigth./ sum(obs_weigth);
            mean_vec = [mean(unit(:,1) .* obs_weigth) mean(unit(:,2) .* obs_weigth)];
            fr_cy = atan2(mean_vec(1), mean_vec(2))/(2*pi);
            
            
        end
        
        function plotSphPatchGrid(el_grid, az_grid, data)
            % Plot a 3D spherical cap 
            % ad hoc function for plotting congruent cells
            % see. Generating Fuhrmann et al. 2013 statistically robust multipath stacking maps using congruent cells
            %
            % used for generateMultipath
            %
            % INPUT
            %   el_grid     array of elevation       [1 x n]
            %   az_grid     cell of array of azimuth [1 x n] of [1 x m(e)]            
            %   data        cell of data elevation x elevation in azimuth [1 x n] of [1 x m(e)]
            %   
            % SYNTAX: 
            %   plotSphPatchGrid(el_grid, az_grid, data)
            
            figure;
            hold on
            el_grid = flipud(el_grid);
            az_grid = fliplr(az_grid);
            data = fliplr(data);
            d_el  = diff(el_grid);
            el_grid = el_grid - [d_el; d_el(end)]/2;
            el_grid =  [el_grid; el_grid(end)+d_el(end)];
            for e = 1 : (length(el_grid)-1)
                data{e} = fliplr(data{e});
                d_az  = diff(az_grid{e}');
                if ~isempty(d_az)
                    az_gridw = az_grid{e}'- [d_az; d_az(end)]/2;
                    az_gridw =  [az_gridw; az_gridw(end)+d_az(end)];
                    for a = 1 : (length(az_gridw)-1)
                        plotSphPatch(el_grid(e:e+1)',az_gridw(a:a+1)', data{e}(a));
                    end
                end
            end
            
            function plotSphPatch(lats, lons, color)
                if lons(2) < lons(1)
                    lons(2) = lons(2)+2*pi;
                end
                dlons = (1 + (1 -cos(mean(lats)))*20)/180*pi;
                lonst = (lons(1):dlons:lons(2))';
                if lonst(end) ~= lons(2)
                    lonst = [lonst; lons(2)];
                end
                lons = lonst;
                lats = pi/2 -lats;
                xyz =[ [cos(lons)*sin(lats(1)) sin(lons)*sin(lats(1)) ones(size(lons))*cos(lats(1))];
                    flipud([cos(lons)*sin(lats(2)) sin(lons)*sin(lats(2)) ones(size(lons))*cos(lats(2))])];
                patch(xyz(:,1),xyz(:,2),xyz(:,3),color);
            end
        end
        
        function fr_cy = circularModedRobustMean(cycle)
            % estimate a roubust mean mean for data over 0 -1
            %
            % SYNTAX:
            %     fr_cy = Core_Utils.circularModedRobustMean(cycle)
            
            % elimnate the nan
            cycle(isnan(cycle)) = [];
            mode_cy = mode(cycle);
            
            % center around the mode and then take the string(robust) mean
            idx_inf = (cycle - mode_cy) < -0.5;
            idx_sup = (cycle - mode_cy) > 0.5;
            cycle(idx_inf) = cycle(idx_inf) +0.5;
            cycle(idx_sup) = cycle(idx_sup) -0.5;
            
            fr_cy = strongMean(cycle,1,0.95,2.5);
            
        end
        
        function amb_idx = getAmbIdx(cycle_slip , obs)
            % get matrix of same dimesion of the observation showing the ambiguity index of the obsarvation
            %
            % SYNTAX:
            % this.getAmbIdx()
            
            amb_idx = ones(size(cycle_slip), 'uint32');
            n_epochs = size(amb_idx,1);
            n_stream = size(amb_idx,2);
            for s = 1:n_stream
                if s > 1
                    amb_idx(:, s) = amb_idx(:, s) + amb_idx(n_epochs, s-1);
                end
                cs = find(cycle_slip(:, s) > 0)';
                for c = cs
                    amb_idx(c:end, s) = amb_idx(c:end, s) + 1;
                end
            end
            amb_idx = amb_idx .* uint32(obs ~= 0);
            amb_idx = Core_Utils.remEmptyAmbIdx(amb_idx);
        end
        
        function createEmptyProject(prj_dir, prj_name, prj_type, prj_site)
            % create empty config file
            %
            % SYNTAX
            %    createEmptyProject(prj_dir, prj_name)
            %    createEmptyProject(prj_name)
            
            fnp = File_Name_Processor();
            
            if nargin == 3
                state = Prj_Settings('', fnp.checkPath(prj_dir));
            else
                state = Prj_Settings('');
            end
            
            if nargin == 1
                prj_name = prj_dir;
                prj_dir = fnp.getFullDirPath([state.getHomeDir filesep '..']);
            end

            if nargin < 4
                prj_site = prj_name;
            end
            
            if iscell(prj_dir)
                
            end
            
            log = Core.getLogger();
            log.addMarkedMessage(sprintf('Creating a new project "%s" into %s', prj_name, prj_dir));
            
            [status, msg, msgID] = mkdir(fnp.checkPath([prj_dir]));
            [status, msg, msgID] = mkdir(fnp.checkPath([prj_dir filesep 'config']));
            [status, msg, msgID] = mkdir(fnp.checkPath([prj_dir filesep 'out']));
            [status, msg, msgID] = mkdir(fnp.checkPath([prj_dir filesep 'out/log']));
            [status, msg, msgID] = mkdir(fnp.checkPath([prj_dir filesep 'RINEX']));
            [status, msg, msgID] = mkdir(fnp.checkPath([prj_dir filesep 'station']));
            [status, msg, msgID] = mkdir(fnp.checkPath([prj_dir filesep 'station/CRD']));
            [status, msg, msgID] = mkdir(fnp.checkPath([prj_dir filesep 'station/ocean']));
            [status, msg, msgID] = mkdir(fnp.checkPath([prj_dir filesep 'station/MET']));
            [status, msg, msgID] = mkdir(fnp.checkPath([prj_dir filesep 'antenna']));
            [status, msg, msgID] = mkdir(fnp.checkPath([prj_dir filesep 'antenna/custom']));
            state.setPrjHome(fnp.checkPath(prj_dir));
            state.prj_name = prj_name;
            state.setOutDir('./out');
            state.setObsDir('./RINEX');
            state.setMetDir('./station/MET');
            state.setCrdDir('station/CRD');
                            
            crd_file = state.getCrdFile;
            if ~exist(crd_file, 'file')
                try
                    % try to create an empty ocean loading file
                    fid = fopen(crd_file, 'w');
                    fwrite(fid, '# Insert here the coordinates of the receivers');
                    fclose(fid);
                catch ex
                end
            end
            
            state.setOceanDir(fullfile(state.getHomeDir, 'station/ocean'));
            state.setOceanFile([prj_site '.blq']);
            ocean_file = state.getOceanFile();
            if ~exist(ocean_file, 'file')
                try
                    % try to create an empty ocean loading file
                    fid = fopen(ocean_file, 'w');
                    fwrite(fid, '$$ Insert here the ocean loading displacements');
                    fclose(fid);
                catch ex
                end
            end

            % Set resources dir:
            data_dir = fnp.getFullDirPath(fullfile(prj_dir, '..', '..'));
            state.eph_dir = fullfile(data_dir, 'satellite', 'EPH', '${WWWW}');
            state.clk_dir = fullfile(data_dir, 'satellite', 'CLK', '${WWWW}');
            state.erp_dir = fullfile(data_dir, 'reference', 'ERP', '${YYYY}');
            state.bias_dir = fullfile(data_dir, 'satellite', 'BIAS');
            state.iono_dir = fullfile(data_dir, 'reference', 'IONO', '${YYYY}');
            state.vmf_dir = fullfile(data_dir, 'reference', 'VMF', '${VMFR}', '${VMFS}', '${YYYY}');            
            state.atm_load_dir = fullfile(data_dir, 'reference', 'ATM_LOAD', '${YYYY}');            
            
            config_path = fnp.checkPath([prj_dir filesep 'config' filesep 'config.ini']);
            state.setSessionDuration(86400);
            state.setBuffer(7200, 7200);
            if nargin >= 3
                switch prj_type
                    case 1 % PPP for tropo
                        state.setToTropoPPP();
                    case 2 % NET no iono, no tropo
                        state.setToShortNET();
                    case 3 % NET no iono
                        state.setToMediumNET();
                        state.rate_ztd_net = 3600;
                    case 4 % NET iono-free
                        state.setToLongNET();
                end
            end
            state.check();
            state.save(config_path);
            Core.getCurrentSettings.import(state);
        end
        
        function y = fillNan1D(y,x)
            % fill the nan into the y
            %
            % SYNTAX
            % y = Core_Utils.fillNan1D(y,<x>)
            
            if nargin < x
                x = 1: length(y);
            end
            idx_nan = isnan(y);
            if sum(~idx_nan) > 2
                int_data = interp1(x(~idx_nan),y(~idx_nan),x(idx_nan));
                y(idx_nan) = int_data;
            end
        end
        
        function [Amp,Phase,f] = getSpectrum(y,smpl_rate)
            % compute the spectrum with fft
            %
            % SYNTAX:
            %  [Amp,Phase,f] = Core_Utils.getSpectrum(y,smpl_rate);
            if nargin < 2
                smpl_rate = 1;
            end
            Y = fft(y);
            
            L = length(y);
            Fs = 1 /smpl_rate;
            f = Fs*(0:(L/2))/L;
            Y = Y(1:(L/2 +1));
            Amp = abs(Y/L)*2;
            Phase = angle(Y);
            
        end
        
        function [val] = spline(t, order)
            % Compute matrix entry for spline
            %
            % INPUT
            %   t -> 0 : 1
            %   order -> 1,3
            %
            % SYNTAX:
            %  Core_Utils.cubicSplic(t)
            switch order
                case 0
                    val = ones(size(t));
                case 1
                    val = Core_Utils.linearSpline(t);
                case 2
                    %%% tBD
                case 3
                    val = Core_Utils.cubicSpline(t);
            end
            
        end
        
        function [id_ko] = snoopGatt(ssat_err, thr, thr_propagate)
            % mark the residual after one thresdol till their mov max  reenter the
            % second threshold
            %
            % SYNTAX:
            %    [w] = Core_Utils.snoopGatt(res, thr, thr_propagate)
            id_ko = false(size(ssat_err));
            for s = 1 : size(id_ko, 2)
                id_ko(:,s) = (movmax(abs(ssat_err(:,s)), 20) > thr_propagate) & flagExpand(abs(ssat_err(:,s)) > thr, 100);
            end
        end
        
        function [idx_ko] = snoopArcLim(ssat_err,thr_propagate)
             % Find begin and end of an arc, if it's out of thr flag it till comes down under threshold 
            %
            % SYNTAX:
            %    [w] = Core_Utils.snoopArcLim(ssat_err,thr_propagate)
            idx_ko = false(size(ssat_err));
            ot = abs(ssat_err) > thr_propagate;
            for s = 1 : size(ssat_err, 2)
                
                % Beginning of the arc
                tmp = ot(:, s) + ~isnan(ssat_err(:, s));
                
                id_start_bad = find(diff([0; tmp]) == 2);
                for i = 1 : numel(id_start_bad)
                    id_stop_bad = find(ot(id_start_bad(i) : end, s) == 0, 1, 'first'); % find when the arc is now under thr
                    
                    % rem over threshold elements
                    if isempty(id_stop_bad)
                        idx_ko(id_start_bad(i) : end, s) = true;
                    else
                        idx_ko(id_start_bad(i) + (0 : (id_stop_bad - 1)), s) = true;
                    end
                end
                
                % End of the arc (flip method)
                ot(:, s) = flipud(ot(:, s));
                tmp = ot(:, s) + flipud(~isnan(ssat_err(:, s)));
                
                id_start_bad = find(diff([0; tmp]) == 2);
                for i = 1 : numel(id_start_bad)
                    id_stop_bad = find(ot(id_start_bad(i) : end, s) == 0, 1, 'first'); % find when the arc is now under thr
                    
                    % rem over threshold elements
                    if isempty(id_stop_bad)
                        idx_ko(size(idx_ko, 1) + 1 - (id_start_bad(i) : size(idx_ko, 1)), s) = true;
                    else
                        idx_ko(size(idx_ko, 1) + 1 - (id_start_bad(i) + (0 : (id_stop_bad - 1))), s) = true;
                    end
                end
            end % ----------------
        end
                
        function [x,inv_diag] = fastInvDiag(N,B,mode)
            % solve the linear system and compute the diagonal entry of the inverse of N square matrix. This is
            % much faster than computing the whole inverse in case of very sparse
            % matrix. Its is done jointly with the resolution of the system in order to
            % keep the decomposition of the matrix.
            % SYNTAX
            %  inv_diag = Core_Utils.fastInvDiag(N,B,<mode>)
            
            if nargin < 3
                mode = 'ldl';
            end
            n_p = size(N,1);
            if strcmpi(mode,'ldl')
                % rememeber : inv(N) ==  P*iL'*inv(D)*iL*P'
                [L,D,P] = ldl(N);
                y = L\B;
                x = (D*L')\y;
                x= P*x;
                iL = inv(L);
                inv_diag = P'*sum(iL.*repmat(1./diag(D),1,n_p).*iL)';
            elseif  strcmpi(mode,'chol')
                [R] = chol(N);
                R = chol(N);
                y = R'\B;
                x=R\y;
                iR = inv(R);
                inv_diag = sum(iR.^2,2);
            end
        end
        
        function [f_nnan] = firstNoNan(mat)
            % find first no nan value of the colums
            %
            % SYNTAX
            %    [f_nnan] = Core_Utils.firstNoNan(mat)
            f_nnan = zeros(1,size(mat,2));
            for i = 1 : length(f_nnan)
                idx = find(~isnan(mat(:,i)),1,'first');
                if ~isempty(idx)
                    f_nnan(i) = mat(idx,i);
                else
                    f_nnan(i) = 0;
                end
            end
        end
        
        function [lid] = ordinal2logical(id,n_el)
            % ordinal 2 logical index
            %
            % SYNTAX
            %    [lid] = Core_Utils.ordinal2logical(id,n_el)
            lid = false(n_el,1);
            lid(id) = true;
        end
        
        
        function x = solveLDLo(L,D,b)
            % solve system where normal matrix has beeen ldl decomposed
            % NOTE : A = L*D*L' (sparse matrices)
            %
            % SYNTAX
            %    x = Core_Utils.solveLDL(L,D,b)
            y = L\b;
            y = y./diag(D);
            x = L'\y;
        end
        
        function x = solveLDL(L,D,b,P,keep_id)
            x = zeros(size(b));
            
            L_red = L(keep_id,keep_id);
            d = diag(D);
            d_red = d(keep_id);
            
            id_est_amb = (P*keep_id)>0;
            
            P_red = P(id_est_amb,keep_id);
            
            x(id_est_amb,:) = P_red * (L_red' \((1./d_red).*(L_red\(P_red' *b(id_est_amb,:)))));
        end
        
        function sys_list = getPrefSys(sys_list)
            % get prefered system (HARDCODED)
            if sum(sys_list == 'G')>0
                sys_list = 'G';
            elseif sum(sys_list == 'R')>0
                sys_list = 'R';
            elseif sum(sys_list == 'E') >0
                sys_list = 'E';
            elseif sum(sys_list == 'C') >0
                sys_list = 'C';
            elseif sum(sys_list == 'J') >0
                sys_list = 'J';
            elseif sum(sys_list == 'I') >0
                sys_list = 'I';
            else
                sys_list = [];
            end
        end
        
        function [dtm, lat, lon, georef, info] = getDTM(nwse, res)
            % Get the dtm of an area delimited by geographical coordinates nwse
            %
            % INPUT
            %   nwse(1)     North [deg: -90:90]
            %   nwse(2)     West  [deg: -180:180]
            %   nwse(3)     South [deg: -90:90]
            %   nwse(4)     East  [deg: -180:180]
            %   res         resolution ('maximum' / 'high' / 'medium' / 'low' / 'v1' / 'v3') - default low
            %
            % For res low 'maximum' to 'low'
            % If the DTM is not found in the DTM folder of the project
            % download it from -> use http://www.marine-geo.org/services/
            %
            % For res low '1' to '3' (they have higher resolutions)
            % If the DTM is not found in the DTM folder 
            % download it from -> use http://viewfinderpanoramas.org/
            %
            % SYNTAX
            %  [dtm, lat, lon, georef, info] = Core_Utils.getDTM(nwse, res)
            
            if nargin <= 1 || isempty(res)
                res = 'low';
            end
                        
            switch res(end)
                case '1'
                    [dtm, lat, lon, georef, info] = Core_Utils.getViewFinderPanoDTM(nwse, '1');
                    % Fallback to lower resolution
                    if ~any(dtm(:))
                        [dtm, lat, lon, georef, info] = Core_Utils.getViewFinderPanoDTM(nwse, '3');
                    end
                case '3'
                    [dtm, lat, lon, georef, info] = Core_Utils.getViewFinderPanoDTM(nwse, '3');
                    % Fallback to lower resolution
                    if ~any(dtm(:))
                        [dtm, lat, lon, georef, info] = Core_Utils.getMarineGeoDTM(nwse, 'high');
                    end
                otherwise
                    [dtm, lat, lon, georef, info] = Core_Utils.getMarineGeoDTM(nwse, res);
            end
        end

        function [dtm, lat, lon, georef, info] = getViewFinderPanoDTM(nwse, res)
            % Get the dtm of an area delimited by geographical coordinates nwse
            %
            % INPUT
            %   nwse(1)     North [deg: -90:90]
            %   nwse(2)     West  [deg: -180:180]
            %   nwse(3)     South [deg: -90:90]
            %   nwse(4)     East  [deg: -180:180]
            %   res         resolution ('1', '3') - high to low
            %
            % If the DTM is not found in the DTM folder of the project
            % download it from -> use http://www.marine-geo.org/services/
            %
            % SYNTAX
            %  [dtm, lat, lon, georef, info] = Core_Utils.getViewFinderPanoDTM(nwse, res)
            
            if nargin < 2 || isempty(res)
                res = '1';
            end
            dtm = zeros(2,2);
            georef = [];
            lat = nwse([3 1]);
            lon = nwse([2 4]);
            info = [];
            %%
            
            if res == '1'
                patch_size = 3601;
            else
                patch_size = 1201;
            end

            % Determine path of storage
            dtm_path = File_Name_Processor.getFullDirPath(fullfile(Core.getInstallDir(), '../data/reference/DTM/VFP1'));            
            % Use local storage dir, because DTM data can be shared among projects
            local_dtm_path = fullfile(Core.getLocalStorageDir, 'reference' , 'DTM', ['VFP' res]);
            if ~(exist(dtm_path, 'dir') == 7)
                dtm_path = local_dtm_path;
            end

            %%
            % Determine the patches to load
            lim_x = sort(floor(nwse([2 4]))); lim_y = sort(floor(nwse([1 3])));
            [X, Y] = meshgrid(lim_x(1):lim_x(2), lim_y(1):lim_y(2));
            X = mod(X+180, 360)-180;
            patch_prefix = 'NSWE';
            patches = [patch_prefix(1 + uint8(Y(:)<0))' num2str(abs(Y(:)), '%02d') patch_prefix(3 + uint8(X(:)>=0))' num2str(abs(X(:)), '%03d')];
            file_list = {};
            folder_list = {};
            file_exist = [];
            for p = 1:size(patches,1)
                folder_list{p} = [iif(Y<0,'S', '') char('A' + (abs(Y(p))/4)) num2str(floor((X(p) + 180)/6) + 1,'%02d')];
                file_list{p,1} = fullfile(dtm_path, folder_list{p}, [patches(p,:) '.hgt']);
                file_exist(p,1) = (exist(file_list{p}, 'file') == 2);
            end

            % Prepare to download missing files
            missing_files = unique(folder_list(~file_exist));
            for i = 1 : numel(missing_files)
                file_to_download = sprintf('viewfinderpanoramas.org/dem%s/%s.zip', res, missing_files{i});
                Core.getLogger.addMarkedMessage(sprintf('Trying to download: "%s",\n this file is big it might require some time depending on your internet connection',file_to_download));
                Core_Utils.downloadHttpTxtResUncompress(file_to_download,dtm_path);
            end
            
            % Recheck file presence
            for p = 1:size(patches,1)
                file_exist(p,1) = (exist(file_list{p}, 'file') == 2);
            end
            
            %%
            if any(file_exist)
                n_y = (diff(minMax(Y))+1);
                n_x = (diff(minMax(X))+1);
                full_dtm = nan(n_y * (patch_size(1)-1) + 1, n_x * (patch_size(end)-1) + 1, 'single');
                lat = linspace(min(Y(:)),max(Y(:))+1, size(full_dtm,1))';
                lon = linspace(min(X(:)),max(X(:))+1, size(full_dtm,2))';
                for p = 1:size(patches,1)
                    % Not considering Earth as a "cylinder"
                    p_x = X(p) - min(X(:)) + 1; % index of the patch in the merged DTM
                    p_y = Y(p) - min(Y(:)) + 1; % index of the patch in the merged DTM
                    file_name = file_list{p};
                    try
                        fid = fopen(file_name,'rb','ieee-be');
                        dtm = single(fread(fid,'*int16'));
                        dtm(dtm == -32768) = nan;
                        dtm = reshape(dtm, patch_size(1), patch_size(end));
                        dtm = rot90(dtm);
                        fclose(fid);
                        full_dtm((p_y-1)*(patch_size(1)-1) + [1:patch_size(1)]', ...
                            (p_x-1)*(patch_size(end)-1) + [1:patch_size(end)]') = dtm;
                    catch ex
                        Core.getLogger.addError(sprintf('Something went wrong reading "%s"', file_name));
                        Core_Utils.printEx(ex)
                    end
                end

                % Cut the DTM as requested
                id_lat = (find(lat < nwse(3), 1, 'last'):find(lat < nwse(1), 1, 'last'))';
                id_lat = max(min(id_lat, numel(lat)), 1);
                id_lon = (find(lon < nwse(2), 1, 'last'):find(lon < nwse(4), 1, 'last'))';
                id_lon = max(min(id_lon, numel(lon)), 1);

                lat = lat(id_lat);
                lon = lon(id_lon);
                dtm = flipud(fillmissing(full_dtm(id_lat, id_lon),'linear','EndValues','none'));
                clear full_dtm;
            end
            
            if ~any(dtm(:))
                Logger.getInstance.addError(sprintf('Failed to download the required DTM file or missing files at this resolution'));
                dtm = zeros(2,2);
                georef = [];
                lat = nwse([3 1]);
                lon = nwse([2 4]);
                info = [];
            end
        end

        function [dtm, lat, lon, georef, info] = getMarineGeoDTM(nwse, res)
            % Get the dtm of an area delimited by geographical coordinates nwse
            %
            % INPUT
            %   nwse(1)     North [deg: -90:90]
            %   nwse(2)     West  [deg: -180:180]
            %   nwse(3)     South [deg: -90:90]
            %   nwse(4)     East  [deg: -180:180]
            %   res         resolution ('maximum' / 'high' / 'medium' / 'low') - default low
            %
            % If the DTM is not found in the DTM folder of the project
            % download it from -> use http://www.marine-geo.org/services/
            %
            % SYNTAX
            %  [dtm, lat, lon, georef, info] = Core_Utils.getMarineGeoDTM(nwse, res)
            
            if nargin < 2 || isempty(res)
                res = 'low';
            end
            
            dtm_path = File_Name_Processor.getFullDirPath(fullfile(Core.getInstallDir(), '../data/reference/DTM/'));
            nwse = [ceil(nwse(1)*1e3)/1e3 floor(nwse(2)*1e3)/1e3 floor(nwse(3)*1e3)/1e3 ceil(nwse(4)*1e3)/1e3]; % Reduce duplicates with small variation of coordinates
            
            local_dtm_path = fullfile(Core.getState.getHomeDir, 'reference' , 'DTM');
            dtm_name = sprintf('dtm_N%06dW%07d_S%06dE%07d_%s.tiff', round(1e2*nwse(1)), round(1e2*nwse(2)), round(1e2*nwse(3)), round(1e2*nwse(4)), res);

            if ~(exist(dtm_path, 'dir') == 7)
                dtm_path = local_dtm_path;
            end
            
            if ~exist(dtm_path, 'dir')
                try
                    mkdir(dtm_path);
                catch ex
                    Logger.getInstance.addWarning(sprintf('Could not create %s\nUsing local folder\nException: %s', dtm_path, ex.message))
                    % no write permission
                    dtm_path = fullfile('reference' , 'DTM');
                    dtm_name = sprintf('dtm_N%06dW%07d_S%06dE%07d_%s.tiff', round(1e2*nwse(1)), round(1e2*nwse(2)), round(1e2*nwse(3)), round(1e2*nwse(4)), res);
                    if ~exist(dtm_path, 'dir')
                        try
                            mkdir(dtm_path);
                        catch ex
                            Logger.getInstance.addError(sprintf('Could not create %s DTM cannot be saved locally\nException: %s\n', dtm_path, ex.message))
                        end
                    end
                end
            end
            [file_info] = dir(fullfile(dtm_path, dtm_name));
            if ~isempty(file_info) && (file_info(1).bytes == 0)
                delete(fullfile(dtm_path, dtm_name)); % remove empty files
            end
            if ~exist(fullfile(dtm_path, dtm_name), 'file')
                if exist(fullfile(local_dtm_path, dtm_name), 'file')
                    dtm_path = local_dtm_path;
                else
                    if ispc()
                        aria2c_path = '.\utility\thirdParty\aria2-extra\aria2_win\aria2c.exe';
                    elseif ismac()
                        aria2c_path = '/usr/local/bin/aria2c';
                        if ~exist(aria2c_path, 'file')
                            aria2c_path = '/opt/homebrew/bin/aria2c';
                        end
                    else % is linux
                        aria2c_path = '/usr/bin/aria2c';
                        if ~exist(aria2c_path, 'file')
                            aria2c_path = '/usr/local/bin/aria2c';
                        end
                    end
                    if exist(aria2c_path, 'file')
                        aria_call = sprintf('%s "%snorth=%f&west=%f&south=%f&east=%f%s%s" --dir="%s" --out="%s"', aria2c_path, 'https://www.gmrt.org/services/GridServer.php?', nwse(1), nwse(2), nwse(3), nwse(4) , '&layer=topo&format=geotiff&resolution=', res, dtm_path, dtm_name);
                        Logger.getInstance.addMarkedMessage(['Executing: "' aria_call '"']);
                        dos(aria_call)
                        [file_info] = dir(fullfile(dtm_path, dtm_name));
                        if ~isempty(file_info) && (file_info(1).bytes == 0)
                            delete(fullfile(dtm_path, dtm_name)); % remove empty files
                        end
                    end
                    if ~exist(aria2c_path, 'file') || ~exist(fullfile(dtm_path, dtm_name), 'file')
                        Logger.getInstance.addWarning('Aria2c failed, try executing websave');
                        url = sprintf('%snorth=%f&west=%f&south=%f&east=%f%s%s', 'https://www.gmrt.org/services/GridServer.php?', nwse(1), nwse(2), nwse(3), nwse(4) , '&layer=topo&format=geotiff&resolution=', res);
                        % Use websave
                        option = weboptions('timeout', 25, 'UserAgent', 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:53.0) Gecko/20100101 Firefox/94.1');
                        res = websave(fullfile(dtm_path, dtm_name), url, option);
                    end

                end
            end
            
            % Read DTM
            try
                [dtm, georef, lat, lon, info] = geotiffReader(fullfile(dtm_path, dtm_name));
            catch ex
                Logger.getInstance.addError(sprintf('Aria failed to download the required DTM file\nException: %s', ex.message));
                dtm = zeros(2,2);
                georef = [];
                lat = nwse([3 1]);
                lon = nwse([2 4]);
                info = [];
            end
        end
        
        function playAlert()
            % Warning at the end of a job! (Play a sound)
            load handel;
            sound(y(1:16000),Fs);
        end
    end
end
