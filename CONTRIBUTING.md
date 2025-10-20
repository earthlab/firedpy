# Contributing to firedpy

Thank you for your interest in contributing to firedpy! This document provides guidelines for both internal Earth Lab team members and external contributors.

## Prerequisites
- Python 3.8 or higher
- Git

## Development Workflows

### For Internal Earth Lab / firedpy Team Members

#### Getting Started
1. **Contact the team**: Reach out to Nate, Aashish, or Adam to request a development branch
2. **Branch creation**: A dev branch will be created for you, cloned from `main`
3. **Clone your dev branch**: `git clone https://github.com/earthlab/firedpy.git -b your-dev-branch`
4. **Set up environment**: Follow the development setup steps below

#### Development Process
1. **Work on your `[your-initials]-dev` branch**: Make all changes on your assigned development branch
2. **Create pull requests**: PR your changes into `nonprod` branch
3. **Testing**: Ensure all tests pass before requesting review
4. **Review process**: Adam or Aashish will review your PRs

#### Branch Naming Convention
- Development branches: `[your-intials]-dev`
- Example: `NH-dev`

#### Branch Workflow Summary
```
[your-initials]-dev (your development branch)
  ↓
nonprod (staging/testing)
  ↓
main (production)
```

### For External Contributors (do we want to do this? better way to set up?)

#### Getting Started
1. **Fork the repository**: Click "Fork" on the GitHub repository page
2. **Clone your fork**: `git clone https://github.com/your-username/firedpy.git`
3. **Set up environment**: Follow the development setup steps below
4. **Create a feature branch**: `git checkout -b feature/your-feature-name`

#### Development Process
1. **Work on your feature branch**: Make changes on your feature branch
2. **Create pull requests**: PR your changes into the `nonprod` branch
3. **Review process**: Internal team members will review your PRs
4. **Testing**: Ensure all tests pass before requesting review


## How to Contribute

### Reporting Issues
Before creating an issue, please:
1. Search existing issues to avoid duplicates
2. Use the an issue template (Bug Report, Feature Request, or Question)
3. Provide as much detail as possible!

### Suggesting Enhancements
When suggesting new features:
1. Use the Feature Request template
2. Provide scientific context and use cases
3. Consider implementation and dependencies
4. Think about validation approaches (if applicable)

### Code Contributions

#### Pull Request Process
1. PR your dev branch / fork into `nonprod` (be sure to a template, either 'major' or 'minor')
2. Tag relevant [issues]
3. Request review from maintainers (Adam or Aashish)

#### Coding Standards
- Follow PEP 8 style guidelines
- Use meaningful variable and function names
- Annotate your code!
- Write tests for new functionality

### Documentation
- Update docstrings for any new functions
- Add usage examples for new features
- Update README if needed

## Scientific Guidelines

### Data Validation
- Test changes + document
- Document any assumptions or limitations

### Reproducibility
- Use version control for all changes
- Document parameter choices
- Provide example workflows
- Consider backward compatibility

## Issue and PR Templates

This repository uses GitHub templates to standardize issue reporting and pull requests:

### Issue Templates
- **Bug Report**: For reporting bugs with scientific context
- **Feature Request**: For suggesting new functionality
- **Question**: For usage help and methodology questions

### Pull Request Templates
- **General**: For most changes (concise template)
- **Major Change/Release**: For major features, breaking changes, or new releases
- **Minor Fix/Update**: For small bug fixes and documentation updates

### Template Selection Guidelines
- **Internal team members**: Use appropriate template based on change scope
- **External contributors**: Use General template for most changes, Major template for significant features

## Getting Help

- Check the documentation first
- Search existing issues
- Use GitHub Discussions for general questions
- Use the Question template for specific help requests

## Recognition

Contributors will be acknowledged in the project documentation and release notes. Significant contributions may warrant co-authorship on relevant publications.

Thank you for contributing to firedpy!
