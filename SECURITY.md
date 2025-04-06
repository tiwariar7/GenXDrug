# Security Policy

## Overview

This document outlines the security measures and policies implemented in the GenXDrug application to protect user data, API keys, and ensure secure operation of the application.

## API Key Security

### Storage
- API keys are stored in a separate `config.py` file
- The `config.py` file is excluded from version control via `.gitignore`
- API keys are never exposed in the application's frontend or logs
- API keys are never transmitted to the client-side

### Usage
- API keys are validated before use
- API keys are cleaned of any whitespace or special characters
- API keys are only used for their intended API endpoints
- API keys are never logged or displayed in error messages

## Data Security

### Input Validation
- All user inputs are validated before processing
- Protein sequences are validated for valid amino acids
- SMILES strings are validated for chemical structure
- File uploads are validated for correct format and size

### Data Transmission
- All API communications use HTTPS
- Sensitive data is never transmitted in plain text
- API responses are validated before processing
- Error messages are sanitized to prevent information leakage

## Application Security

### Dependencies
- All dependencies are regularly updated
- Dependencies are sourced from trusted repositories
- Version numbers are explicitly specified in requirements.txt
- Security vulnerabilities in dependencies are monitored

### Session Management
- Streamlit session state is used for temporary data storage
- Session data is cleared when no longer needed
- No sensitive data is stored in session state
- Session timeouts are implemented where appropriate

## Error Handling

### Logging
- Error messages are sanitized before display
- No sensitive information is included in logs
- Stack traces are only shown in development mode
- API errors are handled gracefully with user-friendly messages

### Monitoring
- API usage is monitored for unusual patterns
- Failed authentication attempts are logged
- Rate limiting is implemented for API calls
- System resources are monitored for potential abuse

## Reporting Security Issues

If you discover a security vulnerability within this application, please follow these steps:

1. Do not disclose the vulnerability publicly until it has been addressed
2. Contact the maintainers immediately at [tiwariar279@gmail.com]
3. Provide detailed information about the vulnerability
4. Include steps to reproduce the issue
5. Wait for acknowledgment and resolution

## Security Updates

- Security patches are applied as soon as they are available
- Users are notified of critical security updates
- Major version updates include security improvements
- Security-related changes are documented in the changelog

## Best Practices

### For Users
- Keep your local environment secure
- Do not share API keys or sensitive data
- Use strong passwords for any authentication
- Report any suspicious activity immediately

### For Developers
- Follow secure coding practices
- Regularly update dependencies
- Review and test security measures
- Document security-related changes

## Compliance

This application follows industry-standard security practices and complies with:
- OWASP security guidelines
- API security best practices
- Data protection regulations
- Secure coding standards

## Contact

For security-related inquiries or to report vulnerabilities, please contact:
- Email: [tiwariar279@gmail.com]
- Security Team: [Team-Propensity]

## Version History

- v1.0.0: Initial security policy implementation
- v1.1.0: Added API key security measures
- v1.2.0: Enhanced error handling and logging
- v1.3.0: Implemented comprehensive input validation 
